#include <fstream>
#include <random>
#include "recon.h"
#include "util.h"
#include "forward.h"
#include "jsonutil.h"
#include "logger.h"
#include "span.h"
#include "amino.h"

Reconstructor::Reconstructor()
  : profileSamples (DefaultProfileSamples),
    profileNodeLimit (0),
    includeBestTraceInProfile (true),
    rndSeed (ForwardMatrix::random_engine::default_seed),
    maxDistanceFromGuide (DefaultMaxDistanceFromGuide),
    usePosteriorsForProfile (true),
    minPostProb (DefaultProfilePostProb)
{ }

bool Reconstructor::parseReconArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-tree") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      treeFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-seqs") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      Require (guideFilename.size() == 0, "Can't specify both -guide and -seqs");
      seqsFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-guide") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      Require (seqsFilename.size() == 0, "Can't specify both -seqs and -guide");
      guideFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-savetree") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      treeSaveFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-saveseqs") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      seqsSaveFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-saveguide") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      guideSaveFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-band") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      maxDistanceFromGuide = atoi (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-noband") {
      maxDistanceFromGuide = -1;
      argvec.pop_front();
      return true;

    } else if (arg == "-samples") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      profileSamples = atoi (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      usePosteriorsForProfile = false;
      return true;

    } else if (arg == "-minpost") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      minPostProb = atof (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      usePosteriorsForProfile = true;
      return true;
      
    } else if (arg == "-states") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      profileNodeLimit = atoi (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-nobest") {
      includeBestTraceInProfile = false;
      argvec.pop_front();
      return true;

    } else if (arg == "-seed") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      rndSeed = atoi (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }

  return diagEnvParams.parseDiagEnvParams (argvec)
    || parseModelArgs (argvec);
}

bool Reconstructor::parseModelArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-model") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      modelFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-savemodel") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      modelSaveFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }

  return false;
}

void Reconstructor::loadReconFiles() {
  Require (seqsFilename.size() > 0 || guideFilename.size() > 0, "Must specify sequences");

  generator = ForwardMatrix::newRNG();
  generator.seed (rndSeed);

  if (modelFilename.size()) {
    LogThisAt(1,"Loading model from " << modelFilename << endl);
    ifstream modelFile (modelFilename);
    ParsedJson pj (modelFile);
    model.read (pj.value);
  } else {
    LogThisAt(1,"Using default amino acid model" << endl);
    model = defaultAminoModel();
  }

  if (seqsFilename.size()) {
    LogThisAt(1,"Loading sequences from " << seqsFilename << endl);
    seqs = readFastSeqs (seqsFilename.c_str());
    if (maxDistanceFromGuide < 0 && treeFilename.size())
      LogThisAt(1,"Don't need guide alignment: banding is turned off and tree is supplied" << endl);
    else {
      LogThisAt(1,"Building guide alignment" << endl);
      AlignGraph ag (seqs, model, 1, diagEnvParams, generator);
      Alignment align = ag.mstAlign();
      guide = align.path;
      gapped = align.gapped();
    }

  } else {
    LogThisAt(1,"Loading guide alignment from " << guideFilename << endl);
    gapped = readFastSeqs (guideFilename.c_str());
    const Alignment align (gapped);
    guide = align.path;
    seqs = align.ungapped;
  }

  if (treeFilename.size()) {
    LogThisAt(1,"Loading tree from " << treeFilename << endl);
    ifstream treeFile (treeFilename);
    tree.parse (JsonUtil::readStringFromStream (treeFile));
  } else {
    LogThisAt(1,"Building neighbor-joining tree" << endl);
    auto dist = model.distanceMatrix (gapped);
    tree.buildByNeighborJoining (gapped, dist);
  }

  if (modelSaveFilename.size()) {
    ofstream modelFile (modelSaveFilename);
    model.write (modelFile);
  }

  if (seqsSaveFilename.size()) {
    ofstream seqsFile (seqsSaveFilename);
    writeFastaSeqs (seqsFile, seqs);
  }

  if (guideSaveFilename.size()) {
    if (gapped.empty())
      Warn("No guide alignment to save");
    else {
      ofstream guideFile (guideSaveFilename);
      writeFastaSeqs (guideFile, gapped);
    }
  }

  if (treeSaveFilename.size()) {
    ofstream treeFile (treeSaveFilename);
    treeFile << tree.toString() << endl;
  }

  buildIndices();
}

void Reconstructor::buildIndices() {
  generator.seed (rndSeed);  // re-seed generator, in case it was used during prealignment

  for (size_t n = 0; n < seqs.size(); ++n) {
    Assert (seqIndex.find (seqs[n].name) == seqIndex.end(), "Duplicate sequence name %s", seqs[n].name.c_str());
    seqIndex[seqs[n].name] = n;
  }

  for (TreeNodeIndex node = 0; node < tree.nodes(); ++node)
    if (!tree.isLeaf(node))
      Assert (tree.nChildren(node) == 2, "Tree is not binary: node %d has %s\nSubtree rooted at %d: %s\n", node, plural(tree.nChildren(node),"child","children").c_str(), node, tree.toString(node).c_str());

  AlignPath reorderedGuide;
  for (TreeNodeIndex node = 0; node < tree.nodes(); ++node)
    if (tree.isLeaf(node)) {
      Assert (tree.nodeName(node).length() > 0, "Leaf node %d is unnamed", node);
      Assert (seqIndex.find (tree.nodeName(node)) != seqIndex.end(), "Can't find sequence for leaf node %s", tree.nodeName(node).c_str());
      const size_t seqidx = seqIndex[tree.nodeName(node)];
      nodeToSeqIndex[node] = seqidx;

      if (!guide.empty())
	reorderedGuide[node] = guide.at(seqidx);

      closestLeaf.push_back (node);
      closestLeafDistance.push_back (0);

      rowName.push_back (tree.seqName(node));

    } else {
      int cl = -1;
      double dcl = 0;
      for (size_t nc = 0; nc < tree.nChildren(node); ++nc) {
	const TreeNodeIndex c = tree.getChild(node,nc);
	const double dc = closestLeafDistance[c] + tree.branchLength(c);
	if (nc == 0 || dc < dcl) {
	  cl = closestLeaf[c];
	  dcl = dc;
	}
      }
      closestLeaf.push_back (cl);
      closestLeafDistance.push_back (dcl);

      rowName.push_back (tree.seqName(node));
    }

  swap (guide, reorderedGuide);
}

void Reconstructor::reconstruct() {
  LogThisAt(1,"Starting reconstruction on " << tree.nodes() << "-node tree" << endl);
  gsl_vector* rootProb = model.insProb;

  LogProb lpFinalFwd = -numeric_limits<double>::infinity(), lpFinalTrace = -numeric_limits<double>::infinity();
  AlignPath path;
  map<int,Profile> prof;
  for (TreeNodeIndex node = 0; node < tree.nodes(); ++node) {
    if (tree.isLeaf(node))
      prof[node] = Profile (model.alphabet, seqs[nodeToSeqIndex[node]], node);
    else {
      const int lChildNode = tree.getChild(node,0);
      const int rChildNode = tree.getChild(node,1);
      const Profile& lProf = prof[lChildNode];
      const Profile& rProf = prof[rChildNode];
      ProbModel lProbs (model, tree.branchLength(lChildNode));
      ProbModel rProbs (model, tree.branchLength(rChildNode));
      PairHMM hmm (lProbs, rProbs, rootProb);

      LogThisAt(2,"Aligning " << lProf.name << " and " << rProf.name << endl);

      ForwardMatrix forward (lProf, rProf, hmm, node, guide.empty() ? GuideAlignmentEnvelope() : GuideAlignmentEnvelope (guide, closestLeaf[lChildNode], closestLeaf[rChildNode], maxDistanceFromGuide));

      Profile& nodeProf = prof[node];
      if (node == tree.root()) {
	path = forward.bestAlignPath();
	nodeProf = forward.bestProfile();
      } else if (usePosteriorsForProfile) {
	BackwardMatrix backward (forward, minPostProb);
	nodeProf = backward.buildProfile (profileNodeLimit, ForwardMatrix::CollapseChains);
      } else
	nodeProf = forward.sampleProfile (generator, profileSamples, profileNodeLimit, ForwardMatrix::CollapseChains, includeBestTraceInProfile);

      const LogProb lpTrace = nodeProf.calcSumPathAbsorbProbs (log_gsl_vector(rootProb), NULL);
      LogThisAt(3,"Forward log-likelihood is " << forward.lpEnd << ", profile log-likelihood is " << lpTrace << " with " << nodeProf.size() << " states" << endl);
      
      if (node == tree.root()) {
	lpFinalFwd = forward.lpEnd;
	lpFinalTrace = lpTrace;
      }

      LogThisAt(5,prof[node].toJson());
    }
  }

  LogThisAt(1,"Final Forward log-likelihood is " << lpFinalFwd << ", final alignment log-likelihood is " << lpFinalTrace << endl);

  vguard<FastSeq> ungapped (tree.nodes());
  for (TreeNodeIndex node = 0; node < tree.nodes(); ++node) {
    if (tree.isLeaf(node)) 
      ungapped[node] = seqs[seqIndex.at(rowName[node])];
    else {
      ungapped[node].seq = string (alignPathResiduesInRow(path.at(node)), '*');
      ungapped[node].name = rowName[node];
    }
  }

  reconstruction = Alignment (ungapped, path);
}

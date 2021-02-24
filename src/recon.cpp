#include <fstream>
#include <random>
#include "recon.h"
#include "util.h"
#include "forward.h"
#include "jsonutil.h"
#include "logger.h"
#include "span.h"
#include "presets.h"
#include "nexus.h"
#include "stockholm.h"
#include "seqgraph.h"
#include "regexmacros.h"
#include "ctok.h"
#include "refiner.h"
#include "memsize.h"
#include "simulator.h"
#include "gamma.h"

const regex nonwhite_re (RE_DOT_STAR RE_NONWHITE_CHAR_CLASS RE_DOT_STAR, regex_constants::basic);
const regex stockholm_re (RE_WHITE_OR_EMPTY "#" RE_WHITE_OR_EMPTY "STOCKHOLM" RE_DOT_STAR);
const regex nexus_re (RE_WHITE_OR_EMPTY "#" RE_WHITE_OR_EMPTY "NEXUS" RE_DOT_STAR);
const regex fasta_re (RE_WHITE_OR_EMPTY ">" RE_DOT_STAR);
//const regex newick_re (RE_WHITE_OR_EMPTY "\\x28" RE_DOT_STAR);  // '('
//const regex json_re (RE_WHITE_OR_EMPTY "\\x7b" RE_DOT_STAR);    // '{'
const regex newick_re (RE_WHITE_OR_EMPTY "\\(" RE_DOT_STAR);  // '('
const regex json_re (RE_WHITE_OR_EMPTY "\\{" RE_DOT_STAR);    // '{'

const vguard<string> Reconstructor::fastAliasArgs = ReconFastAliasArgs;
const vguard<string> Reconstructor::carefulAliasArgs = ReconCarefulAliasArgs;

Reconstructor::Reconstructor()
  : profileSamples (DefaultProfileSamples),
    profileMinLen (0),
    profileMaxLen (numeric_limits<size_t>::max()),
    profileNodeLimit (0),
    maxDPMemoryFraction (DefaultMaxDPMemoryFraction),
    rndSeed (ForwardMatrix::random_engine::default_seed),
    maxDistanceFromGuide (DefaultMaxDistanceFromGuide),
    tokenizeCodons (false),
    guideAlignTryAllPairs (false),
    useUPGMA (true),
    jukesCantorDistanceMatrix (false),
    includeBestTraceInProfile (true),
    keepGapsOpen (false),
    usePosteriorsForProfile (false),
    reconstructRoot (true),
    refineReconstruction (false),
    accumulateSubstCounts (false),
    accumulateIndelCounts (false),
    predictAncestralSequence (false),
    reportAncestralSequenceProbability (false),
    gotPrior (false),
    useLaplacePseudocounts (true),
    usePosteriorsForDot (false),
    useSeparateSubPosteriorsForDot (false),
    keepDotGapsOpen (false),
    minPostProb (0),
    maxEMIterations (DefaultMaxEMIterations),
    minEMImprovement (DefaultMinEMImprovement),
    runMCMC (false),
    outputTraceMCMC (false),
    fixGuideMCMC (false),
    fixTreeMCMC (false),
    fixAlignMCMC (false),
    mcmcSamplesPerSeq (DefaultMCMCSamplesPerSeq),
    mcmcTraceFiles (0),
    outputFormat (StockholmFormat),
    outputLeavesOnly (false),
    guideFile (NULL),
    simulatorRootSeqLen (-1),
    gammaCategories (0),
    gammaShape (1),
    normalizeModel (false)
{ }

int Reconstructor::maxProfileStates() const {
  return profileNodeLimit ? (int) profileNodeLimit : (int) sqrt (maxDPMemoryFraction * getMemorySize() / DPMatrix::cellSize());
}

bool Reconstructor::parseAncSeqArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-ancseq") {
      predictAncestralSequence = true;
      argvec.pop_front();
      return true;

    } else if (arg == "-ancprob") {
      reportAncestralSequenceProbability = true;
      predictAncestralSequence = true;
      argvec.pop_front();
      return true;
    }
  }

  return false;
}

bool Reconstructor::parseReconArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-mcmc") {
      runMCMC = true;
      argvec.pop_front();
      return true;

    } else if (arg == "-savedot") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      dotSaveFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-dotpost") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      minDotPostProb = atof (argvec[1].c_str());
      usePosteriorsForDot = true;
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-dotgapsopen") {
      keepDotGapsOpen = true;
      argvec.pop_front();
      return true;

    } else if (arg == "-dotsubpost") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      minDotSubPostProb = atof (argvec[1].c_str());
      useSeparateSubPosteriorsForDot = true;
      argvec.pop_front();
      argvec.pop_front();
      return true;

    }
  }
  return false;
}

void Reconstructor::checkUniqueSeqFile() {
  Require (fastaGuideFilenames.size() + seqFilenames.size() + nexusGuideFilenames.size() + stockholmGuideFilenames.size() == 1, "Please specify exactly one (and only one) of the following: sequence file, guide alignment, or Nexus file.");
}

void Reconstructor::checkUniqueTreeFile() {
  Require (treeFilename.empty() || (nexusGuideFilenames.empty() && nexusReconFilenames.empty() && stockholmGuideFilenames.empty() && stockholmReconFilenames.empty()), "If you have multiple datasets with trees, please encode each tree in its own Stockholm or Nexus file, rather than specifying the tree file separately.");
  Require (treeFilename.empty() || (seqFilenames.size() + fastaGuideFilenames.size() + (fastaReconFilename.empty() ? 0 : 1) == 1), "If you specify a tree file, there can be one and only one sequence file, otherwise matching up trees to sequence files involves too much guesswork for my liking. To avoid complication, I recommend that if you want to analyze multiple datasets, you please use Nexus or Stockholm format to encode the tree and sequence data directly into the same file.");
}

bool Reconstructor::parseSimulatorArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-rootlen") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      simulatorRootSeqLen = atoi (argvec[1].c_str());
      Require (simulatorRootSeqLen >= 0, "%s must be non-negative", arg.c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-tree") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      simulatorTreeFilenames.push_back (argvec[1]);
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }
  return false;
}

bool Reconstructor::parseModelArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-output") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      const string format = toupper (argvec[1]);
      if (format == "NEXUS")
	outputFormat = NexusFormat;
      else if (format == "FASTA")
	outputFormat = FastaFormat;
      else if (format == "STOCKHOLM")
	outputFormat = StockholmFormat;
      else if (format == "JSON")
	outputFormat = JsonFormat;
      else
	Fail ("Unrecognized format: %s", argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-seed") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      rndSeed = atoi (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-model") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      setModelFilename (argvec[1]);
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-insrate" || arg == "-delrate" || arg == "-insextprob" || arg == "-delextprob"
	       || arg == "-subscale" || arg == "-indelscale" || arg == "-scale"
	       || arg == "-inslen" || arg == "-dellen"
	       || arg == "-gaprate" || arg == "-gaplen" || arg == "-gapextprob") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      const string param = arg.substr(1);
      Require (!modelParam.count(param), "Multiple value for %s specified", arg.c_str());
      modelParam[param] = atof (argvec[1].c_str());
      Require (modelParam[param] >= 0, "%s must be nonnegative", arg.c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-gamma") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      gammaCategories = atoi (argvec[1].c_str());
      if (gammaCategories < 2)
	Warn ("%s must be set to 2 or greater to create discretized-gamma rate categories", arg.c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-shape") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      gammaShape = atof (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-normalize") {
      normalizeModel = true;
      argvec.pop_front();
      return true;

    } else if (arg == "-savemodel") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      modelSaveFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-preset") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      setPresetModelName (argvec[1]);
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-codon") {
      tokenizeCodons = true;
      argvec.pop_front();
      return true;
    }
  }
  return false;
}

bool Reconstructor::parseProfileArgs (deque<string>& argvec, bool allowReconstructions) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-auto") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      const string filename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      switch (detectFormat (filename)) {
      case FastaFormat:
	seqFilenames.push_back (filename);
	break;
      case GappedFastaFormat:
	fastaGuideFilenames.push_back (filename);
	break;
      case NexusFormat:
	if (allowReconstructions) {
	  ifstream in (filename);
	  NexusData nex (in);
	  if (nex.tree.seqNamesBijective(nex.gapped))
	    nexusReconFilenames.push_back (filename);
	  else
	    nexusGuideFilenames.push_back (filename);
	} else
	  nexusGuideFilenames.push_back (filename);
	break;
      case StockholmFormat:
	if (allowReconstructions) {
	  ifstream in (filename);
	  Stockholm stock (in);
	  if (stock.hasTree() && stock.getTree().seqNamesBijective(stock.gapped))
	    stockholmReconFilenames.push_back (filename);
	  else
	    stockholmGuideFilenames.push_back (filename);
	} else
	  stockholmGuideFilenames.push_back (filename);
	break;
      case NewickFormat:
	setTreeFilename (filename);
	break;
      case JsonFormat:
	setModelFilename (filename);
	break;
      case UnknownFormat:
      default:
	Fail ("Could not detect format of file %s; please specify it explicitly", filename.c_str());
	break;
      }
      return true;

    } else if (arg == "-seqs") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      seqFilenames.push_back (argvec[1]);
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-guide") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      fastaGuideFilenames.push_back (argvec[1]);
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-nexus") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      nexusGuideFilenames.push_back (argvec[1]);
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-stockholm") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      stockholmGuideFilenames.push_back (argvec[1]);
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-saveguide") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      guideSaveFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-noancs") {
      outputLeavesOnly = true;
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

    } else if (arg == "-profsamples") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      profileSamples = atoi (argvec[1].c_str());
      usePosteriorsForProfile = false;
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-profminlen") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      profileMinLen = atoi (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-profmaxlen") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      profileMaxLen = atoi (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-profminpost") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      minPostProb = atof (argvec[1].c_str());
      usePosteriorsForProfile = true;
      argvec.pop_front();
      argvec.pop_front();
      return true;
      
    } else if (arg == "-profmaxstates") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      profileNodeLimit = atoi (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-profmaxmem") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      maxDPMemoryFraction = atof (argvec[1].c_str()) / 100;
      profileNodeLimit = 0;
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-nobest") {
      includeBestTraceInProfile = false;
      argvec.pop_front();
      return true;

    } else if (arg == "-keepgapsopen") {
      keepGapsOpen = true;
      argvec.pop_front();
      return true;

    } else if (arg == "-fast") {
      argvec.pop_front();
      for (auto fastArgIter = fastAliasArgs.rbegin(); fastArgIter != fastAliasArgs.rend(); ++fastArgIter)
	argvec.push_front (*fastArgIter);
      return true;

    } else if (arg == "-careful") {
      argvec.pop_front();
      for (auto carefulArgIter = carefulAliasArgs.rbegin(); carefulArgIter != carefulAliasArgs.rend(); ++carefulArgIter)
	argvec.push_front (*carefulArgIter);
      return true;

    } else if (arg == "-rndspan") {
      guideAlignTryAllPairs = false;
      argvec.pop_front();
      return true;

    } else if (arg == "-allspan") {
      guideAlignTryAllPairs = true;
      argvec.pop_front();
      return true;

    } else if (arg == "-upgma") {
      useUPGMA = true;
      argvec.pop_front();
      return true;

    } else if (arg == "-nj") {
      useUPGMA = false;
      argvec.pop_front();
      return true;

    } else if (arg == "-jc") {
      jukesCantorDistanceMatrix = true;
      argvec.pop_front();
      return true;

    } else if (arg == "-tree") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      setTreeFilename (argvec[1]);
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-reroot") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      treeRoot = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }

  return diagEnvParams.parseDiagEnvParams (argvec);
}

bool Reconstructor::parseFitArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-maxiter") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      maxEMIterations = atoi (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-mininc") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      minEMImprovement = atof (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-fixgaprates") {
      accumulateIndelCounts = false;
      argvec.pop_front();
      return true;

    } else if (arg == "-fixsubrates") {
      accumulateSubstCounts = false;
      argvec.pop_front();
      return true;
    }
  }

  return false;
}

bool Reconstructor::parseSamplerArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-samples") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      mcmcSamplesPerSeq = atoi (argvec[1].c_str());
      runMCMC = true;
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-fixguide") {
      fixGuideMCMC = true;
      runMCMC = true;
      argvec.pop_front();
      return true;

    } else if (arg == "-fixtree") {
      fixTreeMCMC = true;
      runMCMC = true;
      argvec.pop_front();
      return true;

    } else if (arg == "-fixalign") {
      fixAlignMCMC = true;
      runMCMC = true;
      argvec.pop_front();
      return true;

    } else if (arg == "-trace") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      mcmcTraceFilename = argvec[1].c_str();
      outputTraceMCMC = true;
      runMCMC = true;
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-norefine") {
      refineReconstruction = false;
      argvec.pop_front();
      return true;

    } else if (arg == "-refine") {
      refineReconstruction = true;
      argvec.pop_front();
      return true;
    }
  }

  return false;
}

bool Reconstructor::parsePremadeArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-recon") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      fastaReconFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-nexusrecon") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      nexusReconFilenames.push_back (argvec[1]);
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-stockrecon") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      stockholmReconFilenames.push_back (argvec[1]);
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }

  return false;
}

bool Reconstructor::parseCountArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-nolaplace") {
      useLaplacePseudocounts = false;
      argvec.pop_front();
      return true;
    } else if (arg == "-refine" || arg == "-norefine") {
      argvec.pop_front();
      return true;
    }
  }

  return false;
}

void Reconstructor::setTreeFilename (const string& fn) {
  Require (treeFilename.empty(), "To specify multiple trees, please encode each one in its own Nexus file, together with the associated sequence data.");
  treeFilename = fn;
}

void Reconstructor::setModelFilename (const string& fn) {
  Require (modelFilename.empty() && presetModelName.empty(), "Please specify one model only.");
  modelFilename = fn;
}

void Reconstructor::setPresetModelName (const string& n) {
  Require (modelFilename.empty() && presetModelName.empty(), "Please specify one model only.");
  presetModelName = n;
}

bool Reconstructor::parseSumArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-counts") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      countFilenames.push_back (argvec[1]);
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }

  return false;
}

void Reconstructor::setModelParam (double& p, const char* paramName) const {
  const string param (paramName);
  if (modelParam.count(param))
    p = modelParam.at(param);
}

void Reconstructor::setModelExtProbByExpectedLength (double& modelProb, const char* lenParamName) const {
  const string param (lenParamName);
  if (modelParam.count(param))
    modelProb = max (0., 1. - 1. / modelParam.at(param));
}

double Reconstructor::getModelParam (const char* paramName, double defaultValue) const {
  const string param (paramName ? paramName : "");
  return (paramName && modelParam.count(param)) ? modelParam.at(param) : defaultValue;
}

void Reconstructor::scaleModel (RateModel& model, const char* subScaleParamName, const char* indelScaleParamName) const {
  const double i = getModelParam(subScaleParamName,1), s = getModelParam(indelScaleParamName,1);
  if (i != 1 || s != 1)
    model = model.scaleRates (i, s);
}

void Reconstructor::loadModel() {
  if (presetModelName.size()) {
    LogThisAt(1,"Loading preset model " << presetModelName << endl);
    model = namedModel (presetModelName);
  } else if (modelFilename.size()) {
    LogThisAt(1,"Loading model from " << modelFilename << endl);
    ifstream modelFile (modelFilename);
    ParsedJson pj (modelFile);
    model.read (pj.value);
  } else if (tokenizeCodons) {
    LogThisAt(1,"Using default codon model (" << DefaultCodonModel << ")" << endl);
    model = namedModel (string (DefaultCodonModel));
  } else {
    LogThisAt(1,"Using default amino acid model (" << DefaultAminoModel << ")" << endl);
    model = namedModel (string (DefaultAminoModel));
  }

  if (normalizeModel)
    model = model.normalizeSubstitutionRate();

  setModelParam (model.insRate, "insrate");
  setModelParam (model.delRate, "delrate");
  setModelParam (model.insExtProb, "insextprob");
  setModelParam (model.delExtProb, "delextprob");
  setModelExtProbByExpectedLength (model.insExtProb, "inslen");
  setModelExtProbByExpectedLength (model.delExtProb, "dellen");

  setModelParam (model.insRate, "gaprate");
  setModelParam (model.delRate, "gaprate");
  setModelParam (model.insExtProb, "gapextprob");
  setModelParam (model.delExtProb, "gapextprob");
  setModelExtProbByExpectedLength (model.insExtProb, "gaplen");
  setModelExtProbByExpectedLength (model.delExtProb, "gaplen");

  scaleModel (model, "subscale", NULL);
  scaleModel (model, NULL, "indelscale");
  scaleModel (model, "scale", "scale");

  if (gammaCategories > 1)
    model = makeDiscretizedGammaModel (model, gammaCategories, gammaShape);

  LogThisAt(2,"Alphabet: " << model.alphabet << endl
	    << "Substitution model has " << plural(model.components(),"mixture component")
	    << ", expected rate " << model.expectedSubstitutionRate() << endl
	    << "Insertion rate " << model.insRate
	    << ", expected insertion length " << model.expectedInsertionLength() << endl
	    << "Deletion rate " << model.delRate
	    << ", expected deletion length " << model.expectedDeletionLength() << endl);

  if (tokenizeCodons)
    codonTokenizer.assertAlphabetTokenized (model.alphabet);

  dataCounts = EventCounts (model, model.components());

  if (modelSaveFilename.size()) {
    ofstream modelFile (modelSaveFilename);
    model.write (modelFile);
  }
}

void Reconstructor::loadTree (Dataset& dataset) {
  Require (treeFilename.size() > 0, "Must specify a tree");
  LogThisAt(1,"Loading tree from " << treeFilename << endl);
  ifstream treeFile (treeFilename);
  dataset.tree.parse (JsonUtil::readStringFromStream (treeFile));
  if (treeRoot.size()) {
    LogThisAt(1,"Re-rooting tree above node " << treeRoot << endl);
    dataset.tree = dataset.tree.rerootAbove (treeRoot);
  }
}

void Reconstructor::buildTree (Dataset& dataset) {
  if (runMCMC && !fixTreeMCMC && !useUPGMA) {
    LogThisAt(1,"Switching to UPGMA tree-estimation algorithm to ensure ultrametric tree for MCMC sampler" << endl);
    useUPGMA = true;
  }
  LogThisAt(1,"Estimating initial tree by " << (useUPGMA ? "UPGMA" : "neighbor-joining") << " (" << dataset.name << ")" << endl);
  auto dist = model.distanceMatrix (dataset.gappedGuide, jukesCantorDistanceMatrix ? 0 : DefaultDistanceMatrixIterations);
  if (useUPGMA)
    dataset.tree.buildByUPGMA (dataset.gappedGuide, dist);
  else
    dataset.tree.buildByNeighborJoining (dataset.gappedGuide, dist);
}

void Reconstructor::seedGenerator() {
  generator = ForwardMatrix::newRNG();
  generator.seed (rndSeed);
}

void Reconstructor::loadSeqs() {
  if (guideSaveFilename.size())
    guideFile = new ofstream (guideSaveFilename);
  for (const auto& fn : seqFilenames)
    loadSeqs (fn, string(), string(), string());
  for (const auto& fn : fastaGuideFilenames)
    loadSeqs (string(), fn, string(), string());
  for (const auto& fn : nexusGuideFilenames)
    loadSeqs (string(), string(), fn, string());
  for (const auto& fn : stockholmGuideFilenames)
    loadSeqs (string(), string(), string(), fn);
  if (guideFile)
    delete guideFile;
}

void Reconstructor::loadSeqs (const string& seqFilename, const string& guideFilename, const string& nexusFilename, const string& stockholmFilename) {
  Require (seqFilename.size() || guideFilename.size() || nexusFilename.size() || stockholmFilename.size(), "Must specify sequences");
  checkUniqueTreeFile();

  if (stockholmFilename.size()) {
    LogThisAt(1,"Loading guide alignment(s) from " << stockholmFilename << endl);
    ifstream stockIn (stockholmFilename);
    while (stockIn && !stockIn.eof()) {
      Stockholm stock (stockIn);
      if (stock.rows() == 0)
	break;
      Dataset& dataset = newDataset();
      dataset.name = stockholmFilename;
      dataset.initGuide (tokenizeCodons ? codonTokenizer.tokenize(stock.gapped) : stock.gapped);
      if (stock.hasTree())
	dataset.tree = stock.getTree();
      else
	buildTree (dataset);
      dataset.prepareRecon (*this);
    }
    
  } else {
    Dataset& dataset = newDataset();

    if (nexusFilename.size()) {
      dataset.name = nexusFilename;
      LogThisAt(1,"Loading guide alignment and tree from " << nexusFilename << endl);
      ifstream nexIn (nexusFilename);
      NexusData nex (nexIn);
      nex.convertNexusToAlignment();
      dataset.tree = nex.tree;
      dataset.initGuide (tokenizeCodons ? codonTokenizer.tokenize(nex.gapped) : nex.gapped);
      dataset.prepareRecon (*this);

    } else {
      if (seqFilename.size()) {
	dataset.name = seqFilename;
	LogThisAt(1,"Loading sequences from " << seqFilename << endl);
	dataset.seqs = readFastSeqs (seqFilename.c_str());
	if (tokenizeCodons)
	  dataset.seqs = codonTokenizer.tokenize (dataset.seqs);
	if (maxDistanceFromGuide < 0 && treeFilename.size())
	  LogThisAt(1,"Don't need guide alignment: banding is turned off and tree is supplied" << endl);
	else {
	  LogThisAt(1,"Building guide alignment (" << dataset.name << ")" << endl);
	  AlignGraph* ag = NULL;
	  if (guideAlignTryAllPairs)
	    ag = new AlignGraph (dataset.seqs, model, 1, diagEnvParams);
	  else {
	    seedGenerator();
	    ag = new AlignGraph (dataset.seqs, model, 1, diagEnvParams, generator);
	  }
	  Alignment align = ag->mstAlign();
	  delete ag;
	  dataset.guide = align.path;
	  dataset.gappedGuide = align.gapped();
	}

      } else {
	LogThisAt(1,"Loading guide alignment from " << guideFilename << endl);
	dataset.name = guideFilename;
	const auto guide = readFastSeqs (guideFilename.c_str());
	dataset.initGuide (tokenizeCodons ? codonTokenizer.tokenize(guide) : guide);
      }

      if (treeFilename.size())
	loadTree (dataset);
      else
	buildTree (dataset);

      dataset.prepareRecon (*this);
    }
  }
}

void Reconstructor::Dataset::initGuide (const vguard<FastSeq>& gapped) {
  gappedGuide = gapped;
  const Alignment align (gappedGuide);
  guide = align.path;
  seqs = align.ungapped;
}

Reconstructor::Dataset& Reconstructor::newDataset() {
  datasets.push_back (Dataset());
  datasets.back().name = string("#") + to_string(datasets.size());
  return datasets.back();
}

void Reconstructor::Dataset::clearPrep() {
  seqIndex.clear();
  nodeToSeqIndex.clear();
  rowName.clear();

  guide.clear();
  closestLeaf.clear();
  closestLeafDistance.clear();
}

void Reconstructor::Dataset::prepareRecon (Reconstructor& recon) {
  tree.validateBranchLengths();

  for (size_t n = 0; n < seqs.size(); ++n) {
    Assert (seqIndex.find (seqs[n].name) == seqIndex.end(), "Duplicate sequence name %s", seqs[n].name.c_str());
    seqIndex[seqs[n].name] = n;
  }

  tree.assertBinary();

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

  if (recon.guideFile && !gappedGuide.empty())
    recon.writeTreeAlignment (tree, gappedGuide, name, *recon.guideFile, false, NULL);
}

void Reconstructor::reconstruct (Dataset& dataset) {
  LogThisAt(1,"Starting reconstruction on " << dataset.tree.nodes() << "-node tree" << " (" << dataset.name << ")" << endl);

  if (!usePosteriorsForProfile)
    seedGenerator();  // re-seed generator, in case it was used during prealignment

  vguard<gsl_vector*> rootProb = model.insProb;
  LogProb lpFinalFwd = -numeric_limits<double>::infinity(), lpFinalTrace = -numeric_limits<double>::infinity();
  const ForwardMatrix::ProfilingStrategy strategy =
    (ForwardMatrix::ProfilingStrategy) (ForwardMatrix::CollapseChains
					| (keepGapsOpen ? ForwardMatrix::KeepGapsOpen : ForwardMatrix::DontKeepGapsOpen)
					| (accumulateSubstCounts ? ForwardMatrix::CountSubstEvents : ForwardMatrix::DontCountSubstEvents)
					| (accumulateIndelCounts ? ForwardMatrix::CountIndelEvents : ForwardMatrix::DontCountIndelEvents)
					| (includeBestTraceInProfile ? ForwardMatrix::IncludeBestTrace : ForwardMatrix::DontIncludeBestTrace));

  SumProduct* sumProd = NULL;
  if (accumulateSubstCounts)
    sumProd = new SumProduct (model, dataset.tree);

  AlignPath path;
  map<int,Profile> prof;
  for (TreeNodeIndex node = 0; node < dataset.tree.nodes(); ++node) {
    if (dataset.tree.isLeaf(node))
      prof[node] = Profile (model.components(), model.alphabet, dataset.seqs[dataset.nodeToSeqIndex[node]], node);
    else {
      const int lChildNode = dataset.tree.getChild(node,0);
      const int rChildNode = dataset.tree.getChild(node,1);
      const Profile& lProf = prof[lChildNode];
      const Profile& rProf = prof[rChildNode];
      ProbModel lProbs (model, dataset.tree.branchLength(lChildNode));
      ProbModel rProbs (model, dataset.tree.branchLength(rChildNode));
      PairHMM hmm (lProbs, rProbs, rootProb);

      LogThisAt(2,"Aligning node #" << lProf.rootRowIndex << " " << lProf.name << " (" << plural(lProf.state.size(),"state") << ", " << plural(lProf.trans.size(),"transition") << ") and node #" << rProf.rootRowIndex << " " << rProf.name << " (" << plural(rProf.state.size(),"state") << ", " << plural(rProf.trans.size(),"transition") << ") to build profile for node #" << node << endl);

      ForwardMatrix* forward = NULL;
      int maxDist = maxDistanceFromGuide;
      while (true) {
	forward = new ForwardMatrix (lProf, rProf, hmm, node, dataset.guide.empty() ? GuideAlignmentEnvelope() : GuideAlignmentEnvelope (dataset.guide, dataset.closestLeaf[lChildNode], dataset.closestLeaf[rChildNode], maxDist), sumProd);
	if (forward->lpEnd > -numeric_limits<double>::infinity())
	  break;
	if (maxDist < 0) {
	  LogThisAt(1,"Sample x-path: (" << to_string_join(forward->x.examplePathToEnd()) << ")\n" << "Sample y-path: (" << to_string_join(forward->y.examplePathToEnd()) << ")\n" << forward->toString(true));
	  hmm.write (clog);
	  Abort ("Zero forward likelihood even in the absence of guide alignment constraints - this is not good");
	}
	if (maxDist*2 > alignPathColumns(dataset.guide)) {
	  LogThisAt(2,"Zero forward likelihood with guide alignment band " << maxDist << "; removing guide alignment constraint" << endl);
	  maxDist = -1;
	} else if (maxDist == 0) {
	  LogThisAt(2,"Zero forward likelihood; increasing guide alignment band from 0 to 1" << endl);
	  maxDist = 1;
	} else {
	  LogThisAt(2,"Zero forward likelihood; doubling guide alignment band from " << maxDist << " to " << (maxDist*2) << endl);
	  maxDist *= 2;
	}
	delete forward;
	forward = NULL;
      }

      if (reconstructRoot)
	LogThisAt(5,"Best alignment of " << lProf.name << " and " << rProf.name << ":\n" << makeAlignmentString (dataset, forward->bestAlignPath(), node, true));

      BackwardMatrix *backward = NULL;
      if (((accumulateSubstCounts || accumulateIndelCounts || !dotSaveFilename.empty()) && node == dataset.tree.root())
	  || (usePosteriorsForProfile && node != dataset.tree.root())) {
	backward = new BackwardMatrix (*forward);
      }

      Profile& nodeProf = prof[node];
      if (node == dataset.tree.root()) {

	if (!dotSaveFilename.empty()) {
	  LogThisAt(3,"Building sequence graph for root node" << endl);
	  const ForwardMatrix::ProfilingStrategy dotStrategy =
	    (ForwardMatrix::ProfilingStrategy)
	    (ForwardMatrix::IncludeBestTrace
	     | (keepDotGapsOpen ? ForwardMatrix::KeepGapsOpen : ForwardMatrix::DontKeepGapsOpen));
	  Profile dotProf = usePosteriorsForDot
	    ? backward->postProbProfile (minDotPostProb, 0, dotStrategy)
	    : backward->bestProfile (dotStrategy);
	  SeqGraph dotSeqGraph (dotProf, model.alphabet, log_vector(model.cptWeight), log_vector_gsl_vector(rootProb), useSeparateSubPosteriorsForDot ? minDotSubPostProb : (usePosteriorsForDot ? minDotPostProb : minPostProb));
	  ofstream dotFile (dotSaveFilename);
	  dotSeqGraph.simplify().writeDot (dotFile);
	}

	if (reconstructRoot) {
	  path = forward->bestAlignPath();
	  nodeProf = forward->bestProfile();
	}
      } else if (usePosteriorsForProfile)
	nodeProf = backward->postProbProfile (minPostProb, maxProfileStates(), strategy);
      else
	nodeProf = forward->sampleProfile (generator, profileSamples, maxProfileStates(), strategy);

      if ((accumulateSubstCounts || accumulateIndelCounts) && node == dataset.tree.root())
	dataset.eigenCounts = backward->getCounts();
      
      if (backward)
	delete backward;
      
      if (node == dataset.tree.root())
	lpFinalFwd = forward->lpEnd;

      if (nodeProf.size()) {
	const LogProb lpTrace = nodeProf.calcSumPathAbsorbProbs (log_vector(model.cptWeight), log_vector_gsl_vector(rootProb), NULL);
	LogThisAt(3,"Forward log-likelihood is " << forward->lpEnd << ", profile log-likelihood is " << lpTrace << " with " << nodeProf.size() << " states" << endl);

	if (node == dataset.tree.root())
	  lpFinalTrace = lpTrace;

	LogThisAt(7,nodeProf.toJson());
      }

      delete forward;
    }
  }

  LogThisAt(2,"Final Forward log-likelihood is " << lpFinalFwd << (reconstructRoot ? (string(", final alignment log-likelihood is ") + to_string(lpFinalTrace)) : string()) << endl);

  if (reconstructRoot) {
    dataset.reconstruction = makeAlignment (dataset, path, dataset.tree.root());
    dataset.gappedRecon = dataset.reconstruction.gapped();

    if (refineReconstruction)
      refine (dataset);
  }

  if (accumulateSubstCounts)
    dataCounts += dataset.eigenCounts.transform (model);
  else if (accumulateIndelCounts)
    dataCounts.indelCounts += dataset.eigenCounts.indelCounts;

  if (sumProd)
    delete sumProd;
}

void Reconstructor::refine (Dataset& dataset) {
  LogThisAt(1,"Refining parent-child alignments (" << dataset.name << ")" << endl);
  vguard<FastSeq>& gappedRecon = dataset.hasAncestralReconstruction() ? dataset.gappedAncestralRecon : dataset.gappedRecon;
  Refiner::History history;
  history.tree = dataset.tree;
  history.gapped = gappedRecon;
  Refiner refiner (model);
  const Refiner::History refinedHistory = refiner.refine (history);
  dataset.tree = refinedHistory.tree;
  gappedRecon = refinedHistory.gapped;
}

void Reconstructor::refineAll() {
  Require (datasets.size() > 0, "Please supply some data");
  for (auto& ds : datasets)
    refine (ds);
}

void Reconstructor::predictAncestors (Dataset& dataset) {
  if (predictAncestralSequence) {
    LogThisAt(1,"Predicting ancestral sequences (" << dataset.name << ")" << endl);
    AlignColSumProduct colSumProd (model, dataset.tree, dataset.gappedRecon);
    while (!colSumProd.alignmentDone()) {
      colSumProd.fillUp();
      colSumProd.fillDown();
      colSumProd.appendAncestralReconstructedColumn (dataset.gappedAncestralRecon);
      if (reportAncestralSequenceProbability)
	colSumProd.appendAncestralPostProbColumn (dataset.gappedAncestralReconPostProb);
      colSumProd.nextColumn();
    }
  }
}

void Reconstructor::predictAllAncestors() {
  for (auto& ds : datasets)
    predictAncestors (ds);
}

void Reconstructor::writeTreeAlignment (const Tree& tree, const vguard<FastSeq>& gapped, const string& name, ostream& out, bool isReconstruction, const ReconPostProbMap* postProb) const {
  Tree t (tree);
  vguard<FastSeq> g (gapped);
  if (outputLeavesOnly) {
    vguard<FastSeq> gl;
    for (TreeNodeIndex n = 0; n < tree.nodes(); ++n)
      if (tree.isLeaf(n))
	gl.push_back (g[n]);
    swap (g, gl);
  }
  if (tokenizeCodons)
    g = codonTokenizer.detokenize (g);
  g = model.convertWildcards (g);
  if (isReconstruction && (outputFormat == NexusFormat || outputFormat == JsonFormat || outputFormat == StockholmFormat)) {
    if (outputLeavesOnly)
      t.assignInternalNodeNames();
    else
      t.assignInternalNodeNames(g);
  }
  switch (outputFormat) {
  case FastaFormat:
    writeFastaSeqs (out, g);
    break;
  case NexusFormat:
    {
      NexusData nexus (g, t);
      nexus.convertAlignmentToNexus();
      nexus.write (out);
    }
    break;
  case JsonFormat:
    writeJson (t, g, out, postProb);
    break;
  case StockholmFormat:
    {
      Stockholm stock (g, t);
      if (postProb) {
	if (outputLeavesOnly)
	  Warn ("Not showing ancestors, so not showing posterior probabilities of ancestors either");
	else
	  for (auto& row_colcharprob: *postProb)
	    for (auto& col_charprob: row_colcharprob.second)
	      for (auto& char_prob: col_charprob.second)
		stock.gs[AncestralSequencePostProbTag][stock.gapped[row_colcharprob.first].name].push_back (string() + to_string(col_charprob.first + 1) + " " + char_prob.first + " " + to_string(char_prob.second));
      }
      stock.gf[StockholmIDTag].push_back (name);
      stock.gf[StockholmLogProbTag].push_back (to_string (TreeAlignFuncs::logLikelihood (model, t, gapped)));
      stock.write (out, 0);
    }
    break;
  default:
    Fail ("Unknown output format");
    break;
  }
}

void Reconstructor::writeJson (const Tree& tree, const vguard<FastSeq>& gapped, ostream& out, const ReconPostProbMap* postProb) const {
  const auto alignCols = gappedSeqColumns (gapped);
  out << "{\"root\": \"" << tree.node[tree.root()].name << "\"," << endl;
  out << " \"branches\": [";
  for (TreeNodeIndex n = 0; n < tree.nodes(); ++n)
    if (n != tree.root())
      out << (n ? "," : "") << "\n  [\"" << tree.node[tree.parentNode(n)].name << "\",\"" << tree.node[n].name << "\"," << tree.node[n].d << "]";
  out << "],\n";
  out << " \"rowData\": {";
  for (int s = 0; s < gapped.size(); ++s) {
    const TreeNodeIndex n = outputLeavesOnly ? tree.findNode (gapped[s].name) : s;
    if (!(!tree.isLeaf(n) && outputLeavesOnly)) {
      out << (s ? "," : "") << "\n  \"" << gapped[s].name << "\": ";
      if (tree.isLeaf(n) || !postProb || !postProb->count(s))
	out << quoted_escaped(gapped[s].seq);
      else {
	out << "[";
	int cols = 0;
	for (auto& col_charprob: postProb->at(s)) {
	  for (; cols < col_charprob.first; ++cols)
	    out << (cols ? "," : "") << "[]";
	  if (cols)
	    out << ",";
	  out << "[";
	  int chars = 0;
	  for (auto& char_prob: col_charprob.second)
	    out << (chars++ ? "," : "") << "[" << quoted_escaped(string(1,char_prob.first)) << "," << to_string(char_prob.second) << "]";
	  out << "]";
	  ++cols;
	}
	for (; cols < alignCols; ++cols)
	  out << (cols ? "," : "") << "[]";
	out << "]";
      }
    }
  }
  out << endl << "}}" << endl;
}

void Reconstructor::writeRecon (const Dataset& dataset, ostream& out) const {
  writeTreeAlignment (dataset.tree, predictAncestralSequence ? dataset.gappedAncestralRecon : dataset.gappedRecon, dataset.name, out, true, reportAncestralSequenceProbability ? &dataset.gappedAncestralReconPostProb : NULL);
}

void Reconstructor::writeRecon (ostream& out) const {
  Assert (datasets.size() > 0, "No dataset");
  for (auto& ds : datasets)
    writeRecon (ds, out);
}

void Reconstructor::writeCounts (ostream& out) const {
  dataCounts.writeJson (out);
}

void Reconstructor::writeModel (ostream& out) const {
  model.write (out);
}

void Reconstructor::loadRecon() {
  if (fastaReconFilename.size()) {

    Dataset& dataset = newDataset();
    dataset.name = fastaReconFilename;

    loadTree (dataset);

    LogThisAt(1,"Loading reconstruction from " << fastaReconFilename << endl);
    dataset.gappedRecon = readFastSeqs (fastaReconFilename.c_str());

    dataset.tree.reorderSeqs (dataset.gappedRecon);
    dataset.reconstruction = Alignment (dataset.gappedRecon);

    dataset.gappedGuide = dataset.gappedRecon;
  }

  for (const auto& nexusReconFilename : nexusReconFilenames) {
    Dataset& dataset = newDataset();
    dataset.name = nexusReconFilename;

    LogThisAt(1,"Loading reconstruction and tree from " << nexusReconFilename << endl);

    ifstream nexIn (nexusReconFilename);
    NexusData nex (nexIn);
    nex.convertNexusToAlignment();
    dataset.tree = nex.tree;
    dataset.gappedRecon = nex.gapped;

    dataset.tree.reorderSeqs (dataset.gappedRecon);
    dataset.reconstruction = Alignment (dataset.gappedRecon);

    dataset.gappedGuide = dataset.gappedRecon;
  }

  for (const auto& stockholmReconFilename : stockholmReconFilenames) {

    LogThisAt(1,"Loading reconstructions and trees from " << stockholmReconFilename << endl);

    ifstream stockIn (stockholmReconFilename);
    size_t nStock = 0;
    while (stockIn && !stockIn.eof()) {
      Stockholm stock (stockIn);
      if (stock.rows() == 0)
	break;
      Require (stock.hasTree(), "Stockholm alignment lacks tree");
      Dataset& dataset = newDataset();
      dataset.name = stockholmReconFilename + " alignment #" + to_string(++nStock);
      dataset.gappedRecon = stock.gapped;
      dataset.tree = stock.getTree();
      dataset.tree.reorderSeqs (dataset.gappedRecon);
      dataset.reconstruction = Alignment (dataset.gappedRecon);
      dataset.gappedGuide = dataset.gappedRecon;
    }
  }
}

void Reconstructor::loadCounts() {
  if (countFilenames.empty())
    priorCounts = EventCounts (model, model.components());
  else
    for (auto iter = countFilenames.begin(); iter != countFilenames.end(); ++iter) {
      ifstream in (*iter);
      ParsedJson pj (in);
      EventCounts c;
      c.read (pj.value);
      if (iter == countFilenames.begin())
	priorCounts = c;
      else
	priorCounts += c;
      gotPrior = true;
    }
  if (useLaplacePseudocounts) {
    priorCounts += EventCounts (priorCounts, priorCounts.components(), 1.);
    gotPrior = true;
  }
  dataCounts = priorCounts;
}

void Reconstructor::count (Dataset& dataset) {
  dataset.eigenCounts = EigenCounts (model.components(), model.alphabetSize());
  dataset.eigenCounts.accumulateCounts (model, dataset.reconstruction, dataset.tree, accumulateIndelCounts, accumulateSubstCounts);
  if (accumulateSubstCounts)
    dataCounts += dataset.eigenCounts.transform (model);
  else if (accumulateIndelCounts)
    dataCounts.indelCounts += dataset.eigenCounts.indelCounts;
}

Reconstructor::HistoryLogger::HistoryLogger (Reconstructor& recon, const string& name)
  : recon (&recon),
    out (NULL),
    name (name)
{
  if (recon.outputTraceMCMC && recon.mcmcTraceFilename.size())
    out = new ofstream (recon.mcmcTraceFilename + "." + to_string(++recon.mcmcTraceFiles));
}

Reconstructor::HistoryLogger::~HistoryLogger() {
  if (out)
    delete out;
}

void Reconstructor::HistoryLogger::logHistory (const Sampler::History& history) {
  if (recon->outputTraceMCMC)
    recon->writeTreeAlignment (history.tree, history.gapped, name, out ? *out : cout, true);
}

void Reconstructor::sampleAll() {
  Require (datasets.size() > 0, "Please supply some data");
  Require (!fixAlignMCMC || !fixTreeMCMC, "You can't fix both tree and alignment when doing MCMC - you must sample one of them!");
  if (runMCMC) {
    SimpleTreePrior treePrior;
    vguard<Sampler> samplers;
    vguard<HistoryLogger*> loggers;
    size_t totalNodes = 0;
    CachingRateModel cachedModel (model);
    for (auto& dataset: datasets) {
      if (!dataset.hasReconstruction())
	reconstruct (dataset);
      if (!dataset.hasAncestralReconstruction())
	predictAncestors (dataset);
      vguard<FastSeq>& gappedRecon = dataset.hasAncestralReconstruction() ? dataset.gappedAncestralRecon : dataset.gappedRecon;
      dataset.tree.assignInternalNodeNames (gappedRecon);
      samplers.push_back (Sampler (cachedModel, treePrior, dataset.gappedGuide));
      loggers.push_back (new HistoryLogger (*this, dataset.name));
      Sampler& sampler = samplers.back();
      sampler.addLogger (*loggers.back());
      sampler.useFixedGuide = fixGuideMCMC;
      sampler.sampleAncestralSeqs = dataset.hasAncestralReconstruction();
      Sampler::History history;
      history.tree = dataset.tree;
      history.gapped = gappedRecon;
      sampler.initialize (history, dataset.name);
      if (fixTreeMCMC)
	sampler.fixTree();
      if (fixAlignMCMC)
	sampler.fixAlignment();
      totalNodes += history.tree.nodes();
    }

    const unsigned int nSamples = mcmcSamplesPerSeq * totalNodes;
    LogThisAt(1,"Starting MCMC sampler ("
	      << plural(mcmcSamplesPerSeq,"sample") << " per node, "
	      << plural(nSamples,"sample") << " in total)" << endl);
    Sampler::run (samplers, generator, nSamples);

    for (size_t n = 0; n < datasets.size(); ++n) {
      Dataset& dataset = datasets[n];
      Sampler& sampler = samplers[n];
      dataset.tree = sampler.bestHistory.tree;
      dataset.gappedRecon = sampler.bestHistory.gapped;
      dataset.reconstruction = Alignment (dataset.gappedRecon);
      dataset.clearPrep();

      if (refineReconstruction)
	refine (dataset);
    }

    for (HistoryLogger* logger: loggers)
      delete logger;
  }
}

void Reconstructor::reconstructAll() {
  Require (datasets.size() > 0, "Please supply some data");
  for (auto& ds : datasets)
    reconstruct (ds);
}

void Reconstructor::countAll() {
  Require (datasets.size() > 0, "Please supply some data");
  dataCounts = EventCounts (model, model.components());
  for (auto& ds : datasets)
    if (ds.hasReconstruction())
      count (ds);
    else
      reconstruct (ds);
  dataPlusPriorCounts = dataCounts + priorCounts;
}

void Reconstructor::fit() {
  Require (accumulateIndelCounts || accumulateSubstCounts, "With indel AND substitution rates fixed, model has no free parameters to fit.");
  if (datasets.empty()) {
    Require (gotPrior, "Please specify some data, or pseudocounts, in order to fit a model.");
    priorCounts.optimize (model, accumulateIndelCounts, accumulateSubstCounts);
  } else {
    LogProb lpLast = -numeric_limits<double>::infinity();

    priorCounts.indelCounts.lp = 0;
    for (size_t iter = 0; iter < maxEMIterations; ++iter) {
      countAll();
      const LogProb lpData = dataCounts.indelCounts.lp, lpPrior = gotPrior ? priorCounts.logPrior (model, accumulateIndelCounts, accumulateSubstCounts) : 0;
      const LogProb lpWithPrior = lpData + lpPrior;
      LogThisAt (1, "EM iteration #" << iter + 1 << ": log-likelihood" << (gotPrior ? (string(" (") + to_string(lpData) + ") + log-prior (" + to_string(lpPrior) + ")") : string()) << " = " << lpWithPrior << endl);
      if (lpWithPrior <= lpLast + abs(lpLast)*minEMImprovement)
	break;
      const LogProb oldExpectedLogLike = dataCounts.expectedLogLikelihood (model) + lpPrior;
      dataPlusPriorCounts.optimize (model, accumulateIndelCounts, accumulateSubstCounts);
      const LogProb newExpectedLogLike = dataCounts.expectedLogLikelihood (model) + (gotPrior ? priorCounts.logPrior (model, accumulateIndelCounts, accumulateSubstCounts) : 0);
      LogThisAt(5, "Expected log-likelihood went from " << oldExpectedLogLike << " to " << newExpectedLogLike << " during M-step" << endl);
      lpLast = lpWithPrior;
    }
  }
}

Alignment Reconstructor::makeAlignment (const Dataset& dataset, const AlignPath& path, TreeNodeIndex root) const {
  vguard<FastSeq> ungapped (dataset.tree.nodes());
  for (TreeNodeIndex node : dataset.tree.nodeAndDescendants(root)) {
    if (dataset.tree.isLeaf(node))
      ungapped[node] = dataset.seqs[dataset.seqIndex.at(dataset.rowName[node])];
    else {
      ungapped[node].seq = string (alignPathResiduesInRow(path.at(node)), Alignment::wildcardChar);
      ungapped[node].name = dataset.rowName[node];
    }
  }
  return Alignment (ungapped, path);
}

string Reconstructor::makeAlignmentString (const Dataset& dataset, const AlignPath& path, TreeNodeIndex root, bool assignInternalNodeNames) const {
  vguard<FastSeq> g = makeAlignment(dataset,path,root).gapped();
  for (TreeNodeIndex node = 0; node < dataset.tree.nodes(); ++node)
    if (g[node].name.empty())
      g[node].name = dataset.tree.seqName(node);
  Tree tbig = dataset.tree;
  if (assignInternalNodeNames)
    tbig.assignInternalNodeNames (g);
  Tree t (tbig.toString (root));
  vguard<FastSeq> gt;
  for (auto n : dataset.tree.nodeAndDescendants(root))
    gt.push_back (g[n]);
  Stockholm stock (gt, t);
  ostringstream out;
  stock.write (out, 0);
  return out.str();
}

Reconstructor::FileFormat Reconstructor::detectFormat (const string& filename) {
  LogThisAt(3,"Auto-detecting format for file " << filename << endl);
  ifstream in (filename);
  if (!in)
    Fail ("File not found: %s", filename.c_str());
  string line;
  do {
    if (in.eof())
      Fail ("Couldn't auto-detect file format (all whitespace): %s", filename.c_str());
    getline(in,line);
  } while (!regex_match (line, nonwhite_re));

  if (regex_match (line, stockholm_re)) {
    LogThisAt(3,"Detected Stockholm format" << endl);
    return StockholmFormat;
  } else if (regex_match (line, nexus_re)) {
    LogThisAt(3,"Detected Nexus format" << endl);
    return NexusFormat;
  } else if (regex_match (line, newick_re)) {
    LogThisAt(3,"Detected Newick format" << endl);
    return NewickFormat;
  } else if (regex_match (line, json_re)) {
    LogThisAt(3,"Detected JSON format" << endl);
    return JsonFormat;
  } else if (regex_match (line, fasta_re)) {
    in.close();
    const vguard<FastSeq> seqs = readFastSeqs (filename.c_str());
    for (auto& fs : seqs)
      for (char c : fs.seq)
	if (Alignment::isGap(c)) {
	  LogThisAt(3,"Detected gapped FASTA format" << endl);
	  return GappedFastaFormat;
	}
    LogThisAt(3,"Detected FASTA format" << endl);
    return FastaFormat;
  }

  LogThisAt(3,"Format unknown" << endl);
  return UnknownFormat;
}

void Reconstructor::simulate() {
  Require (simulatorTreeFilenames.size(), "Please provide a tree");
  if (simulatorRootSeqLen < 0) {
    LogThisAt (1, "Using default root sequence length of " << DefaultSimulatorRootSeqLen << endl);
    simulatorRootSeqLen = DefaultSimulatorRootSeqLen;
  }
  loadModel();
  seedGenerator();
  for (const auto& simulatorTreeFilename: simulatorTreeFilenames) {
    LogThisAt(1,"Loading tree from " << simulatorTreeFilename << endl);
    ifstream treeFile (simulatorTreeFilename);
    if (!treeFile) {
      Warn ("Couldn't open %s", simulatorTreeFilename.c_str());
      continue;
    }
    Tree tree (JsonUtil::readStringFromStream (treeFile));
    if (!tree.nodes()) {
      Warn ("Tree %s is empty", simulatorTreeFilename.c_str());
      continue;
    }
    if (outputFormat != FastaFormat)
      tree.assignInternalNodeNames();
    Stockholm stock = Simulator::simulateTree (generator, model, tree, simulatorRootSeqLen);
    if (tokenizeCodons) {
      stock.gapped = codonTokenizer.detokenize (stock.gapped);
      // for now, just throw away component annotation for codon models
      stock.gr.clear();
    }
    const auto& id = simulatorTreeFilename;
    stock.gf[StockholmIDTag].push_back (id);
    if (outputFormat == StockholmFormat)
      stock.write (cout, 0);
    else
      writeTreeAlignment (tree, stock.gapped, id, cout, false, NULL);
  }
}

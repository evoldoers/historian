#include <gsl/gsl_math.h>
#include "refiner.h"
#include "recon.h"

#define REFINER_EPSILON 1e-3
#define REFINER_NEAR_EQ(X,Y) (gsl_fcmp (X, Y, REFINER_EPSILON) == 0)

#define max3(X,Y,Z) max(max(X,Y),Z)

Refiner::BranchMatrix::BranchMatrix (const RateModel& model, const PosWeightMatrix& xSeq, const PosWeightMatrix& ySeq, TreeBranchLength dist, const GuideAlignmentEnvelope& env, const vguard<SeqIdx>& xEnvPos, const vguard<SeqIdx>& yEnvPos, AlignRowIndex x, AlignRowIndex y)
  : BranchMatrixBase (model, xSeq, ySeq, dist, env, xEnvPos, yEnvPos, x, y)
{
  ProgressLog (plog, 5);
  plog.initProgress ("Branch alignment matrix (%u*%u)", xSize, ySize);

  lpStart() = 0;
  for (SeqIdx xpos = 0; xpos < xSize; ++xpos) {

    plog.logProgress (xpos / (double) (xSize - 1), "row %d/%d", xpos + 1, xSize);

    for (SeqIdx ypos = 0; ypos < ySize; ++ypos)
      if (inEnvelope (xpos, ypos)) {
	XYCell& dest = xyCell (xpos, ypos);

	if (xpos > 0 && inEnvelope (xpos - 1, ypos)) {
	  const XYCell& dSrc = xyCell (xpos - 1, ypos);

	  dest(ProbModel::Delete) = max3 (dSrc(ProbModel::Match) + md,
					  dSrc(ProbModel::Insert) + id,
					  dSrc(ProbModel::Delete) + dd);
	}

      	if (ypos > 0 && inEnvelope (xpos, ypos - 1)) {
	  const XYCell& iSrc = xyCell (xpos, ypos - 1);
	  const LogProb yEmitScore = yEmit[ypos - 1];

	  dest(ProbModel::Insert) = yEmitScore + max (iSrc(ProbModel::Match) + mi,
						      iSrc(ProbModel::Insert) + ii);
	}

	if (xpos > 0 && ypos > 0 && inEnvelope (xpos - 1, ypos - 1)) {
	  const XYCell& mSrc = xyCell (xpos - 1, ypos - 1);
	  const LogProb xyEmitScore = logMatch (xpos, ypos);
	  
	  dest(ProbModel::Match) = xyEmitScore + max3 (mSrc(ProbModel::Match) + mm,
						       mSrc(ProbModel::Insert) + im,
						       mSrc(ProbModel::Delete) + dm);
	}
      }
  }

  const XYCell& endCell = xyCell (xSize - 1, ySize - 1);

  lpEnd = max3 (endCell(ProbModel::Match) + me,
		endCell(ProbModel::Insert) + ie,
		endCell(ProbModel::Delete) + de);

  if (LoggingThisAt(9))
    writeToLog(9);

  LogThisAt(6,"Viterbi log-likelihood is " << lpEnd << endl);
}

AlignPath Refiner::BranchMatrix::best() const {
  CellCoords coords (xSize - 1, ySize - 1, ProbModel::End);
  AlignRowPath xPath, yPath;
  while (coords.xpos > 0 || coords.ypos > 0) {
    bool x, y;
    getColumn (coords, x, y);
    if (x || y) {
      xPath.push_back (x);
      yPath.push_back (y);
    }
    CellCoords src = coords;
    if (x) --src.xpos;
    if (y) --src.ypos;
    const XYCell& srcCell = xyCell (src.xpos, src.ypos);
    const LogProb e = lpEmit (coords);
    LogProb bestSrcLogProb = -numeric_limits<double>::infinity();
    CellCoords bestSrc = src;
    bool foundBest = false;
    for (src.state = 0; src.state < (unsigned int) ProbModel::End; ++src.state) {
      const LogProb srcLogProb = srcCell(src.state) + lpTrans ((State) src.state, (State) coords.state) + e;
      if (srcLogProb > bestSrcLogProb) {
	bestSrcLogProb = srcLogProb;
	bestSrc = src;
	foundBest = true;
      }
    }

    Assert (foundBest, "Could not find traceback state from cell %s", coords.toString().c_str());
    Assert (REFINER_NEAR_EQ (bestSrcLogProb, cell(coords)), "Traceback score (%g) doesn't match stored value (%g) at cell %s", bestSrcLogProb, cell(coords), coords.toString().c_str());

    coords = bestSrc;
  }
  AlignPath path;
  path[xRow] = AlignRowPath (xPath.rbegin(), xPath.rend());
  path[yRow] = AlignRowPath (yPath.rbegin(), yPath.rend());

  LogThisAt(9,"Optimal parent-child alignment:" << endl << alignPathString(path));

  return path;
}

Refiner::Refiner (const RateModel& model)
  : model (model),
    maxDistanceFromGuide (DefaultMaxDistanceFromGuide)
{ }

GuideAlignmentEnvelope Refiner::makeGuide (const Tree& tree, const AlignPath& path, TreeNodeIndex node1, TreeNodeIndex node2) const {
  return GuideAlignmentEnvelope (path, node1, node2, maxDistanceFromGuide);
}

vguard<SeqIdx> Refiner::guideSeqPos (const AlignPath& path, AlignRowIndex row) const {
  return TreeAlignFuncs::getGuideSeqPos (path, row, row);
}

Refiner::History Refiner::refine (const History& oldHistory, TreeNodeIndex node) const {
  const TreeNodeIndex parent = oldHistory.tree.parentNode (node);

  LogThisAt(4,"Attempting branch refinement move between...\n   node #" << node << ": " << oldHistory.tree.seqName(node) << "\n parent #" << parent << ": " << oldHistory.tree.seqName(parent) << endl);

  const TreeBranchLength dist = oldHistory.tree.branchLength(parent,node);
  
  const Alignment oldAlign (oldHistory.gapped);
  const AlignPath oldBranchPath = branchPath (oldAlign.path, oldHistory.tree, node);
  const GuideAlignmentEnvelope newBranchEnv = makeGuide (oldHistory.tree, oldBranchPath, parent, node);

  const AlignPath pCladePath = cladePath (oldAlign.path, oldHistory.tree, parent, node);
  const AlignPath nCladePath = cladePath (oldAlign.path, oldHistory.tree, node, parent);
  
  const vguard<SeqIdx> parentEnvPos = guideSeqPos (oldAlign.path, parent);
  const vguard<SeqIdx> nodeEnvPos = guideSeqPos (oldAlign.path, node);
  
  map<TreeNodeIndex,TreeNodeIndex> exclude;
  exclude[node] = parent;
  exclude[parent] = node;

  const auto pwms = getConditionalPWMs (model, oldHistory.tree, oldHistory.gapped, exclude, allExceptNodeAndAncestors(oldHistory.tree,parent), nodeAndAncestors(oldHistory.tree,parent));
  const PosWeightMatrix& pSeq = pwms.at (parent);
  const PosWeightMatrix& nSeq = pwms.at (node);

  const BranchMatrix branchMatrix (model, pSeq, nSeq, dist, newBranchEnv, parentEnvPos, nodeEnvPos, parent, node);
  const AlignPath newBranchPath = branchMatrix.best();

  const vguard<AlignPath> mergeComponents = { pCladePath, newBranchPath, nCladePath };
  const AlignPath newPath = alignPathMerge (mergeComponents);

#ifdef DEBUG
  LogThisAt(12,"Test of conditional probability weight matrix calculation:" << endl << branchConditionalDump(model, oldHistory.tree, oldHistory.gapped, parent, node));
#endif /* DEBUG */

  LogThisAt(7,"Old (parent:node) alignment:" << endl << alignPathString(oldBranchPath)
	    << "Log-likelihood: " << branchMatrix.logPathProb(oldBranchPath) << endl);

  LogThisAt(6,"New (parent:node) alignment:" << endl << alignPathString(newBranchPath)
	    << "Log-likelihood: " << branchMatrix.logPathProb(newBranchPath) << endl);

  LogThisAt(7,"New full alignment:" << endl << alignPathString(newPath));

  const Alignment newAlign (oldAlign.ungapped, newPath);

  History newHistory;
  newHistory.tree = oldHistory.tree;
  newHistory.gapped = newAlign.gapped();
  
  return newHistory;
}

Refiner::History Refiner::refine (const History& oldHistory) const {
  const Tree& tree = oldHistory.tree;
  tree.assertPostorderSorted();
  History bestHistory = oldHistory;
  LogProb bestLogProb = logLikelihood (bestHistory);
  TreeNodeIndex node = 0;
  int stepsSinceImprovement = 0;
  while (stepsSinceImprovement < tree.nodes() - 1) {
    const History newBestHistory = refine (bestHistory, node);
    const LogProb newBestLogProb = logLikelihood (newBestHistory);
    if (newBestLogProb > bestLogProb) {
      LogThisAt(3,"Branch refinement improved alignment log-likelihood from " << bestLogProb << " to " << newBestLogProb << endl);
      bestHistory = newBestHistory;
      bestLogProb = newBestLogProb;
      stepsSinceImprovement = 0;
    } else {
      if (newBestLogProb < bestLogProb && !REFINER_NEAR_EQ(newBestLogProb,bestLogProb))
	Warn ("During branch refinement, alignment log-likelihood dropped from %g to %g", bestLogProb, newBestLogProb);
      ++stepsSinceImprovement;
      LogThisAt(4,"Branch refinement failed to improve alignment log-likelihood for " << plural(stepsSinceImprovement,"step") << endl);
    }
    node = (node + 1) % (tree.nodes() - 1);  // skip root
  }
  return bestHistory;
}

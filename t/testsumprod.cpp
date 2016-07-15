#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/model.h"
#include "../src/jsonutil.h"
#include "../src/sumprod.h"
#include "../src/logger.h"

int main (int argc, char **argv) {
  if (argc != 4) {
    cout << "Usage: " << argv[0] << " <model> <alignment> <tree>\n";
    exit (EXIT_FAILURE);
  }
  
  RateModel rates;
  ifstream in (argv[1]);
  ParsedJson pj (in);
  rates.read (pj.value);

  vguard<FastSeq> gapped = readFastSeqs (argv[2]);

  ifstream treeStream (argv[3]);
  Tree tree (JsonUtil::readStringFromStream (treeStream));

  tree.reorderSeqs (gapped);
  
  //  logger.setVerbose (8);
  
  AlignColSumProduct sp (rates, tree, gapped);

  while (!sp.alignmentDone()) {
    sp.fillUp();
    sp.fillDown();

    cout << "Column #" << sp.col << endl;
    for (auto node : sp.ungappedRowIndices())
      if (node != sp.columnRoot()) {
	const TreeNodeIndex parent = tree.parentNode(node);
	for (int cpt = 0; cpt < rates.components(); ++cpt)
	  for (AlphTok a = 0; a < rates.alphabetSize(); ++a)
	    for (AlphTok b = 0; b < rates.alphabetSize(); ++b)
	      cout << "P( " << tree.seqName(parent) << " = " << rates.alphabet[a] << cpt << " , " << tree.seqName(node) << " = " << rates.alphabet[b] << cpt << " ) = " << exp (sp.logBranchPostProb (cpt, node, a, b)) << endl;
      }
    vguard<LogProb> lnpp = sp.logNodePostProb (sp.columnRoot());
    for (AlphTok a = 0; a < rates.alphabetSize(); ++a)
      cout << "P( " << tree.seqName(sp.columnRoot()) << " = " << rates.alphabet[a] << " ) = " << exp (lnpp[a]) << endl;
    cout << endl;
    
    sp.nextColumn();
  }

  exit (EXIT_SUCCESS);
}

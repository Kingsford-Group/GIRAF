#include <sstream>
#include <iostream>
#include <fstream>

#include "label_types.h" 
#include "bigraph.h"
#include "scored_set.h"
int non_star = 0;
#include "mica.h"
#include "timer.h"

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 4) {

    cerr << "Usage: " << argv[0] << " <work-directory> <threshold> <type>" << endl << endl;
    cerr << "type = 0 for all bicliques" << endl;
    cerr << "       1 for non-star bicliques" << endl;
    cerr << "       2 for just edges" << endl;
    return 0;
  }

  Labelled_Bigraphs<left_label_t, right_label_t, edge_label_t> graphs;
  double threshold = atof(argv[2]);
  non_star = atoi(argv[3]);

  ofstream results((string(argv[1])+"/bicliques").c_str());

  if(!results) {

    cerr << "Error in opening results file" << endl;
    return 0;
  }

  read_trees((string(argv[1])+"/left_trees").c_str(), num_of_left_trees, left_trees);
  read_trees((string(argv[1])+"/right_trees").c_str(), num_of_right_trees, right_trees);
  graphs.read((string(argv[1])+"/graph.labelled").c_str());

  Timer T;
  cout << T.start() << endl;

  for(set<edge_label_t>::iterator it = graphs.getEdgeTypes().begin();
      it != graphs.getEdgeTypes().end(); it++) {

    ostringstream output;

    if(non_star == 2)
      good_edges<left_label_t, right_label_t, tree_label_t>(graphs.getGraph(*it), threshold, results);
    else
      mica<left_label_t, right_label_t, tree_label_t>(graphs.getGraph(*it), threshold, output, cerr);

    if(!output.str().empty()) {

      results << "Label: " << *it << endl << endl << output.str() << endl;
      results << "------------" << endl;
    }
  }

  cout << T.stop() << endl;
  cout << T.report() << endl;

  results.close();

  return 0;
}

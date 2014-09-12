#ifndef _SCORED_SET_H_
#define _SCORED_SET_H_

#include <map>
#include <set>
#include <fstream>
#include <algorithm>
#include "label_types.h"

int num_of_left_trees; int num_of_right_trees;
map<left_label_t, set<tree_label_t> > left_trees;
map<right_label_t, set<tree_label_t> > right_trees;


template<class node_label_t, class tree_label_t>
void read_trees(const char* filename, int& num_of_trees, map<node_label_t, set<tree_label_t> >& trees) {

  ifstream input_file(filename);

  if(!input_file) {

    cerr << "Could not read the file " << filename << endl;
    input_file.close();
    return;
  }

  string mark; int num_of_lines;
  input_file >> mark >> num_of_trees >> num_of_lines;
  for(int i = 1; i <= num_of_lines; i++) {

    node_label_t node; int num_of_entries; tree_label_t entry;
    input_file >> node >> num_of_entries;
    for(int j = 1; j <= num_of_entries; j++) {

      input_file >> entry; trees[node].insert(entry);
    }
  }

  input_file.close();
}

template<class node_label_t, class tree_label_t>
class scored_set {

 private:

  set<node_label_t> nodes;
  set<tree_label_t> trees;
  
 public:

  scored_set(const set<node_label_t>& given_nodes, const set<tree_label_t>& given_trees)
    : nodes(given_nodes), trees(given_trees) { }

  scored_set(const set<node_label_t>& given_nodes, 
	     map<node_label_t, set<tree_label_t> >& tree_map, int num_of_trees) : nodes(given_nodes) {

    typename set<node_label_t>::iterator it = nodes.begin();
    trees = tree_map[*it];
  
    for(it++; it != nodes.end(); it++) {

      set<node_label_t> result;
      set_union(trees.begin(), trees.end(), 
		tree_map[*it].begin(), tree_map[*it].end(),
		inserter(result, result.begin()));
      swap(result, trees);
    }

    set<node_label_t> temp;
    for(int i = 0; i <= num_of_trees; i++)
      if(trees.find(i) == trees.end())
	temp.insert(i);

    swap(temp, trees);
  }

  const set<node_label_t>& get_nodes() const { return nodes; }
  bool operator<(const scored_set& a) const 
    { return nodes < a.nodes; }

  double score(int num_of_trees) const { 
    return (nodes.size() > 0 ? (num_of_trees-trees.size())/double(num_of_trees) : 0); }

  scored_set intersect(const scored_set& to_intersect) const {

      set<node_label_t> node_intersection;
      set_intersection(nodes.begin(), nodes.end(), to_intersect.nodes.begin(), to_intersect.nodes.end(),
		       inserter(node_intersection, node_intersection.begin()));

      set<tree_label_t> tree_intersection;
      set_union(trees.begin(), trees.end(), to_intersect.trees.begin(), to_intersect.trees.end(),
		inserter(tree_intersection, tree_intersection.begin()));

      return scored_set(node_intersection, tree_intersection);
  }
};

#endif

#ifndef _BIGRAPH_H_
#define _BIGRAPH_H_

#include <iostream>
#include <map>		
#include <set>		
#include <iterator>	
#include <fstream>
#include <sstream>

using namespace std;

template<class left_label_t, class right_label_t>
class Bigraph {

 private:

  map<left_label_t, set<right_label_t> > left_adjacency;
  map<right_label_t, set<left_label_t> > right_adjacency;
  
 public:

  map<left_label_t, set<right_label_t> >& get_left_adjacency() { return left_adjacency; }
  map<right_label_t, set<left_label_t> >& get_right_adjacency() { return right_adjacency; }

  size_t read(const char* filename);
  size_t remove_biclique(const set<left_label_t>& left_nodes, const set<right_label_t>& right_nodes);
};

template<class left_label_t, class right_label_t>
size_t Bigraph<left_label_t, right_label_t>::read(const char* filename) {

  ifstream input_file(filename);

  if (!input_file) {

    cerr << "Could not read the file " << filename << endl;
    input_file.close();
    return 0;
  }

  size_t edge_count(0); left_label_t left_node; right_label_t right_node;
  while(input_file >> left_node >> right_node) {

    left_adjacency[left_node].insert(right_node);
    right_adjacency[right_node].insert(left_node);
    edge_count++;
  }

  input_file.close();
  return edge_count;
}

template<class left_label_t, class right_label_t>
size_t Bigraph<left_label_t, right_label_t>::remove_biclique(const set<left_label_t>& left_nodes, const set<right_label_t>& right_nodes) {

  size_t removed(0);

  for(typename set<left_label_t>::iterator left_it = left_nodes.begin(); left_it != left_nodes.end(); left_it++)
    for(typename set<right_label_t>::iterator right_it = right_nodes.begin(); right_it != right_nodes.end(); right_it++) {
      left_adjacency[*left_it].erase(*right_it);
      right_adjacency[*right_it].erase(*left_it);
    }

  for(typename set<left_label_t>::iterator left_it = left_nodes.begin(); left_it != left_nodes.end(); left_it++)
    if(left_adjacency[*left_it].size() == 0) {
      left_adjacency.erase(*left_it);
      removed++;
    }

  for(typename set<right_label_t>::iterator right_it = right_nodes.begin(); right_it != right_nodes.end(); right_it++)
    if(right_adjacency[*right_it].size() == 0) {
      right_adjacency.erase(*right_it);
      removed++;
    }

  return removed;
}

template<class left_label_t, class right_label_t, class edge_label_t>
class Labelled_Bigraphs {

 private:

  map<edge_label_t, Bigraph<left_label_t, right_label_t> > bigraphs;
  set<edge_label_t> edge_types;

 public:

  void read(const char* filename) {

    ifstream input_file(filename);
    
    if (!input_file) {

      cerr << "Could not read the file " << filename << endl;
      input_file.close();
      return;
    }

    left_label_t left_node; right_label_t right_node;
    char line[1000];
    
    while(input_file >> left_node >> right_node) {

      input_file.getline(line, 1000);
      istringstream line_stream(line);

      edge_label_t edge_type;
      while(line_stream >> edge_type) {
	
	edge_types.insert(edge_type); 
	bigraphs[edge_type].get_left_adjacency()[left_node].insert(right_node);
	bigraphs[edge_type].get_right_adjacency()[right_node].insert(left_node);
      }
    }

    input_file.close();
  }

  Bigraph<left_label_t, right_label_t>& getGraph(edge_label_t edge_type) { return bigraphs[edge_type]; } 
  set<edge_label_t>& getEdgeTypes() { return edge_types; }
 
};

#endif

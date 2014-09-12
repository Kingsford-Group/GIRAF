#ifndef _MICA_H_
#define _MICA_H_

//
// This code is based on an implementation by Wen-Chieh Chang that can be found at
// http://genome.cs.iastate.edu/supertree/download/biclique/
//

#include <iostream>
#include <algorithm>	
#include <fstream>
#include "bigraph.h"
#include "scored_set.h"

const double SIZE = 50;

template<class data_t>
void print_set(const set<data_t>& given_set, ostream& output_stream) {

  typename set<data_t>::iterator it = given_set.begin(); output_stream << *it; 
  for(it++; it != given_set.end(); it++)
    output_stream << " " << *it;
}

template<class in_label_t, class out_label_t>
set<out_label_t> retrieve_other_set(const set<in_label_t>& input_set, map<in_label_t, set<out_label_t> >& adjacency) {

  set<out_label_t> output_set;

  typename set<in_label_t>::iterator it = input_set.begin();
  output_set = adjacency[*it];

  for(it++; it != input_set.end(); it++) {

    set<out_label_t> result;
    set_intersection(output_set.begin(), output_set.end(), 
		     adjacency[*it].begin(), adjacency[*it].end(),
		     inserter(result, result.begin()));
    swap(result, output_set);
  }

  return output_set;
}

template<class left_label_t, class right_label_t, class tree_label_t>
void output_results(ostream& results, Bigraph<left_label_t, right_label_t>& graph, double threshold, const set<scored_set<right_label_t, tree_label_t> >& all) {

  for(typename set<scored_set<right_label_t, tree_label_t> >::iterator it = all.begin(); it != all.end(); it++) {

    scored_set<left_label_t, tree_label_t> 
      left_nodes(retrieve_other_set<right_label_t, left_label_t>(it->get_nodes(), graph.get_right_adjacency()), left_trees, num_of_left_trees);

    if(left_nodes.score(num_of_left_trees) > threshold && 
       (non_star == 0 || (left_nodes.get_nodes().size() > 1 && it->get_nodes().size() > 1))) {

      print_set<left_label_t>(left_nodes.get_nodes(), results); results << endl;
      print_set<right_label_t>(it->get_nodes(), results); results << endl;
      results << left_nodes.score(num_of_left_trees) << " " << it->score(num_of_right_trees) << endl << endl;
    }
  }
}

template<class left_label_t, class right_label_t, class tree_label_t>
void expand_stars(Bigraph<left_label_t, right_label_t>& graph, double threshold, const set<scored_set<right_label_t, tree_label_t> >& stars, 
		  set<scored_set<right_label_t, tree_label_t> >& result) {

  result.clear();

  for(typename set<scored_set<right_label_t, tree_label_t> >::iterator it1 = stars.begin(); it1 != stars.end(); it1++)
    for(typename set<scored_set<right_label_t, tree_label_t> >::iterator it2 = it1; ++it2 != stars.end(); ) {

      scored_set<right_label_t, tree_label_t> intersection = it1->intersect(*it2);

      if(intersection.score(num_of_right_trees) > threshold && stars.find(intersection) == stars.end() &&       
	 result.find(intersection) == result.end())
	  result.insert(intersection);
    }
}

template<class left_label_t, class right_label_t, class tree_label_t>
void expand(Bigraph<left_label_t, right_label_t>& graph, double threshold, const set<scored_set<right_label_t, tree_label_t> >& stars, 
	    const set<scored_set<right_label_t, tree_label_t> >& current, const set<scored_set<right_label_t, tree_label_t> >& all, 
	    set<scored_set<right_label_t, tree_label_t> >& result) {

  result.clear();

  for(typename set<scored_set<right_label_t, tree_label_t> >::iterator it1 = stars.begin(); it1 != stars.end(); it1++)
    for(typename set<scored_set<right_label_t, tree_label_t> >::iterator it2 = current.begin(); it2 != current.end(); it2++) {

      scored_set<right_label_t, tree_label_t> intersection = it1->intersect(*it2);
      if(intersection.score(num_of_right_trees) > threshold && all.find(intersection) == all.end() && 
	 result.find(intersection) == result.end())
	  result.insert(intersection);
    }
}

template<class left_label_t, class right_label_t, class tree_label_t>
size_t mica(Bigraph<left_label_t, right_label_t>& graph, double threshold, ostream& results, ostream & debug = cout) {

  size_t rank(1);
  set<scored_set<right_label_t, tree_label_t> > stars, current, all, result;

  for(typename map<left_label_t, set<right_label_t> >::iterator it = graph.get_left_adjacency().begin(); it != graph.get_left_adjacency().end(); it++) {
    
    scored_set<right_label_t, tree_label_t> star(it->second, right_trees, num_of_right_trees);
    if(star.score(num_of_right_trees) > threshold) {
      stars.insert(star);
    }
  }

  all = stars;

  rank++;
  expand_stars(graph, threshold, stars, current);
  copy(current.begin(), current.end(), inserter(all, all.begin()));

  while (!current.empty()) {

    rank++;
    expand(graph, threshold, stars, current, all, result); swap(current, result);
    copy(current.begin(), current.end(), inserter(all, all.begin()));
  }

  output_results(results, graph, threshold, all);
  
  return all.size();
}

template<class left_label_t, class right_label_t, class tree_label_t>
  size_t good_edges(Bigraph<left_label_t, right_label_t>& graph, double threshold, ostream& results, ostream & debug = cout) {

  size_t count = 0;
  
  for(typename map<left_label_t, set<right_label_t> >::iterator it = graph.get_left_adjacency().begin(); it != graph.get_left_adjacency().end(); it++) {
    
    set<left_label_t> left_node; left_node.insert(it->first);
    scored_set<left_label_t, tree_label_t> 
      scored_left_node(left_node, left_trees, num_of_left_trees);
    
    if(scored_left_node.score(num_of_left_trees) > threshold) {
      
      for(typename set<right_label_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
	
	set<right_label_t> right_node; right_node.insert(*it2);
	scored_set<right_label_t, tree_label_t> 
	  scored_right_node(right_node, right_trees, num_of_right_trees);
	
	if(scored_right_node.score(num_of_right_trees) > threshold) {
	  
	  results << it->first << endl;
	  results << *it2 << endl << endl;
	  results << scored_left_node.score(num_of_left_trees) << " " << scored_right_node.score(num_of_right_trees) << endl << endl;
	  count++;
	}
      }
    }
  }
  
  return count;
}

#endif

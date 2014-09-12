#ifndef TREE_H
#define TREE_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <vector>
#include <set>
#include "util.h"

using namespace std;


//
// Represents a tree. The tree is represented by a pointer to
// its root node.
//
struct TreeNode 
{
  string id;                       // name of node
  double length;                   // length of edge to parent
  int depth;                       // distance from root

  vector<TreeNode *> children;     // children of node
  TreeNode * parent;               // ptr to parent
  vector<TreeNode *> leaves;       // leaves of this node

  // iterators over the node's children
  typedef vector<TreeNode *>::const_iterator const_child_iterator;
  typedef vector<TreeNode *>::iterator child_iterator;

  // Constructor 
  TreeNode(const string & i = "") : id(i), length(0.0),parent(0) {}
};

//
// Read a tree in NH format
//
TreeNode * ReadTree(istream &);

typedef map<string, string> NodeNameMapping;

void ReadNexTranslate(istream &, NodeNameMapping *);

void ReadNexTrees(istream &, NodeNameMapping *, vector<TreeNode *> &, int = 0);

//
// Options that can be passed to WriteTree() to control how
// the tree is output.
//
enum WriteTreeOptions { InternalLabels = 1, NoLengths = 2, InternalIDs = 4, SameLine = 8 } ;

//
// Write a tree in NH format; The output options are
// encoded as an | of the above enum
//
void WriteTree(ostream &, const TreeNode *, int);

// assign valid lists to the leaves member
void AssignLeavesLists(TreeNode *);

// return the names of the leaves as strings
set<string> LeavesOf(const TreeNode *);

// Free the memory associated with a tree
void DeleteTree(TreeNode *);

// Give each non-labeled node an id.
void AssignIDs(TreeNode *, int *);

// make sure we have valid edge lengths
void CheckFixLengths(TreeNode *, double);

void WriteTreeAsDot( ostream &, TreeNode *);

typedef map<string, map<string, double> > DistanceMatrix;

vector<DistanceMatrix> * AllDistanceMatrices( vector<TreeNode *> & trees); 

void ComputeDistanceMatrix(TreeNode *T, DistanceMatrix & M);

void ScaleDistanceMatrix(TreeNode *, DistanceMatrix &);

void PrintDistances(ofstream &, vector<DistanceMatrix> &);


#endif

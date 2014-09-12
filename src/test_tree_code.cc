#include "tree.h"
#include "splits.h"

int
main(int argc, char * argv[])
{
    if (argc < 2)
    {
        cerr << "Usage: test_tree_code trees.nex" << endl << endl;

        cerr << "Applys the first traslate command found in the .nex file to" << endl
             << "every tree record in the file. It then writes out the mapping and" << endl
             << "the trees." << endl;
        exit(3);
    }
    ifstream nexus(argv[1]);

    // read the leaf name mapping if there is one
    NodeNameMapping leafs, *ptr;
    ReadNexTranslate(nexus, &leafs);
    ptr = (leafs.empty()?0:&leafs);
    if (!ptr)
    {
        WARN("No translate table found!");
    }

    // reset the file
    nexus.seekg(0, ios::beg);

    // write what we found
    cout << "Found " << leafs.size() << " mapping entries:" << endl;
    PrintMap(cout, leafs, " -> ", "\n");

    // read the tree collection
    vector<TreeNode *> trees;
    ReadNexTrees(nexus, ptr, trees);

    // write what we found
    cout << "Found " << trees.size() << " trees" << endl;

    for(vector<TreeNode *>::iterator Tree = trees.begin();
        Tree != trees.end();
        ++Tree)
    {
        WriteTree(cout, *Tree, SameLine);
        cout << endl;
    }


    // find all the splits
    SplitDatabase splits;
    AllSplits(trees, splits);

    cout << "Found " << splits.size() << " splits." << endl;

    PrintSplitsReadable(cout, splits);

    PrintSplitsMapping(cout, splits);
}

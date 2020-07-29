#ifndef MASSDAWGNODE_H
#define MASSDAWGNODE_H

#include <vector>
#include <string>
#include <iostream>

using namespace std;

class MassDawgNode{
public:
    // kmer values associated with incoming mass values
    vector<string> kmers;
    // outgoing mass edges
    vector<MassDawgNode *> children;
    // the sinlgy and doubly mass of this node
    float singlyMass;
    float doublyMass;

    // empty constructor
    MassDawgNode ();

    // init with masses and a string
    MassDawgNode (float singlyMass, float doublyMass, string kmer);

    ~MassDawgNode();

    /**
     * Add a kmer to the set of kmers. No duplicates will be made
     * 
     * @param kmer  string      kmer to add to node
    */
    void addKmer (string kmer);

    /**
     * Add a child node to the node called on by creating a connecting edge
     * 
     * @param singlyMass    float  singly charged mass to connect nodes
     * @param doublyMass    float  doubly charged mass to connect nodes
     * @param kmer          string  the string to associate with the new kmer
     * 
     * @return MassDawgNode *   the new child added
    */
    MassDawgNode * addChild(float singlyMass, float doublyMass, string kmer);

    /**
     * Turn the node's value into a large string in order to make 
     * it hashable
     * 
     * @return string   the hashable string
    */
   string hash();

   /**
    * Recursively show this node and all subsequent nodes and edges
    * 
    * @param spaces int     the number of spaces to prepend before printing
   */
    void show(int spaces);
};
#endif
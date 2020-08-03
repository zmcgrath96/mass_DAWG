#ifndef MASSDAWGC_H
#define MASSDAWGC_H

#define MIN(a, b)       (a > b ? b : a)    

#include <unordered_map>
#include <vector>
#include <list>
#include <iostream>
#include <string>

using namespace std;

/**
 * Used for qsort on floats. If d1 < d2, a number < 0 returned
 * if d1 > d2, a number > 0 returned
 * 
*/
int dblCmp(const void * d1, const void * d2);

/**
 * Convert the ppm value of a given mass to a value in daltons
 * 
 * @param mass      float  the mass to calculate the dalton tolerance for
 * @param ppmTol    int     the tolerance in ppm
 * 
 * @return float   the tolerance in daltons 
*/
float ppmToDa(float mass, int ppmTol);

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

class UncheckedNode {
public:
    MassDawgNode * parent;
    MassDawgNode * child;

    UncheckedNode(){}

    ~UncheckedNode(){}
};

class PreviousSequence {
public: 
    vector<float> singlySequence;
    vector<float> doublySequence;
    vector<MassDawgNode *> nodes;

    PreviousSequence() {}

    ~PreviousSequence(){}
};

class MassDawg{
public:
    // empty constructor takes no values
    MassDawg();

    // destructor 
    ~MassDawg();

    /**
     * Show the graph as a tree in the console
    */
   void show();

   /**
     * Add a new singly and doubly charged sequence associated with the kmer to the graph
     * 
     * @param singlySequence    vector<float>  the singly charged sequence of masses
     * @param doublySequence    vector<float>  the doubly charged sequence of masses
     * @param kmer              string          the sequence of amino acids associated with this mass
    */
    void insert(vector<float> singlySequence, vector<float> doublySequence, string kmer);

    /**
     * Search for the input sequence while allowing for up to gapAllowances
     * before the search returns however deep it is in the graph
     * 
     * @param sequence      vector<float>  the sequence to search 
     * @param gapAllowance  int             The number of gaps to allow in the search
     * @param ppmTol        int             the tolerance in parts per million to accept when searching
     * 
     * @return vector<string>               All kmers that we found in the search
    */
   vector<string> fuzzySearch(vector<float> sequence, int gapAllowance, int ppmTol);

    /**
     * Any remaining unchecked nodes will be checked for merging to 
     * complete the dawg. 
    */
    void finish();

private:

    list<UncheckedNode> uncheckedNodes;
    unordered_map<string, MassDawgNode *> minimizedNodes;
    PreviousSequence ps;
    MassDawgNode * root;    

    /**
     * What makes this a graph and not a tree. Combines nodes that share edges and values
     * 
     * @param downTo    int     the level down to which we need to minimize 
    */
    void minimize(int downTo);

    /**
     * Checks to see if the new sequences are greater than the old previous sequence
     * 
     * @param singlySequence    vector<float>  the new sequence of singly charged masses
     * @param doublySequence    vector<float>  the new sequence of doubly charged masses
     * 
     * @return bool     True if the new sequences are greater than the prvious, False otherwise
    */
    bool previousIsLessThan(vector<float> singlySequence, vector<float> doublySequence);

    /**
     * Recursive search of the graph allowing for gapAllowance missed masses in the
     * search before returning whatever is found at the level
     * 
     * @param sequence      vector<float>  The sequence to use to navigate the graph
     * @param currentNode   MassDawgNode *  The current node to investigate
     * @param currentGap    int             The number of gaps we have allowed up until this point
     * @param gapAllowance  int             The total number of gaps to allow
     * @param ppmTol        int             the tolerance in parts per million to accept when searching
     * 
     * @return vector<string>   The kmers associated with the deepest part of the branch investigated
    */
    vector<string> fuzzySearchRec(vector<float> sequence, MassDawgNode * currentNode, int currentGap, int gapAllowance, int ppmTol);

};
#endif
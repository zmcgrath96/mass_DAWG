#ifndef MASSDAWG_H
#define MASSDAWG_H

#define MIN(a, b)       (a > b ? b : a)    

#include <unordered_map>
#include <list>
#include <vector>
#include <iostream>

#include "MassDawgNode.hpp"

using namespace std;

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

class LongestCommonPrefix {
public:
    vector<float> singlySequence;
    vector<float> doublySequence;
    vector<MassDawgNode *> nodes;

    LongestCommonPrefix() {}
    LongestCommonPrefix(vector<float> sS, vector<float> dS, vector<MassDawgNode *> nodes);

    ~LongestCommonPrefix() {}
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

    /**
     * Find the longest common prefix of input sequences to a path in the tree. Used for out
     * of order insertions.
     * 
     * @param singlySequence    vector<float>   the singly sequence we are looking for a prefix of
     * @param doublySequence    vector<float>   the doubly sequence we are looking for a prefix of
     * 
     * @returns LongestCommonPrefix   the class instance holding the longest common prefixe
    */
    LongestCommonPrefix longestCommonPrefix(vector<float> singlySequence, vector<float> doublySequence);
};
#endif
#ifndef MASSDAWG_H
#define MASSDAWG_H

#define MIN(a, b)       (a > b ? b : a)    

#include <unordered_map>
#include <list>
#include <vector>
#include <iostream>

#include "MassDawgNode.hpp"
#include "Edge.hpp"

using namespace std;

class UncheckedNode {
public:
    MassDawgNode * parent;
    MassDawgNode * child;
    Edge * connector;

    UncheckedNode(){}

    ~UncheckedNode(){}
};

class PreviousSequence {
public: 
    vector<double> singlySequence;
    vector<double> doublySequence;
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
     * @param singlySequence    vector<doubly>  the singly charged sequence of masses
     * @param doublySequence    vector<doubly>  the doubly charged sequence of masses
     * @param kmer              string          the sequence of amino acids associated with this mass
    */
    void insert(vector<double> singlySequence, vector<double> doublySequence, string kmer);

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
     * @param singlySequence    vector<double>  the new sequence of singly charged masses
     * @param doublySequence    vector<double>  the new sequence of doubly charged masses
     * 
     * @return bool     True if the new sequences are greater than the prvious, False otherwise
    */
    bool previousIsLessThan(vector<double> singlySequence, vector<double> doublySequence);

};
#endif
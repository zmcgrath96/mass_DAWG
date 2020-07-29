#include "MassDawgNode.hpp"
#include "Edge.hpp"
#include "utils.cpp"


MassDawgNode::MassDawgNode (){}

// init with a string
MassDawgNode::MassDawgNode (string kmer){
        this->kmers.push_back(kmer);
    }

MassDawgNode::~MassDawgNode(){}

/**
 * Add a kmer to the set of kmers. No duplicates will be made
 * 
 * @param kmer  string      kmer to add to node
*/
void MassDawgNode::addKmer (string kmer){
    // check to see if this kmer exists in the set. If not, add it
    for (int i = 0; i < this->kmers.size(); i++){
        if (this->kmers[i].compare(kmer) == 0) return;
    }
    this->kmers.push_back(kmer);
}

/**
 * Add a child node to the node called on by creating a connecting edge
 * 
 * @param singlyMass    double  singly charged mass to connect nodes
 * @param doublyMass    double  doubly charged mass to connect nodes
 * @param kmer          string  the string to associate with the new kmer
 * 
 * @return Edge *       edge connecting the parent to the new child
*/
Edge * MassDawgNode::addChild(double singlyMass, double doublyMass, string kmer){
    Edge * edge = new Edge(singlyMass, doublyMass, kmer);
    this->edges.push_back(edge);

    return edge;
}

/**
 * Turn the node's value into a large string in order to make 
 * it hashable
 * 
 * @return string   the hashable string
*/
string MassDawgNode::hash(){
    string hashable = "";

    // if we dont have any edges, return ""
    if (this->edges.size() == 0) return "";


    // go throught the egdes and add them to a list to sort such that all 
    // nodes will hash the same if they are the same
    double toSortMasses[this->edges.size() * 2];
    for (int i = 0; i < this->edges.size(); i++){
        toSortMasses[i] = this->edges[i]->singlyMass; 
        toSortMasses[i+1] = this->edges[i]->doublyMass;
    }
    qsort(toSortMasses, this->edges.size(), sizeof(double), dblCmp);

    // go through each of these values and add them to hashable
    for (int i = 0; i < this->edges.size(); i ++){
        hashable += to_string(toSortMasses[i]) + to_string(toSortMasses[i+1]);
    }

    return hashable;
}

/**
* Recursively show this node and all subsequent nodes and edges
* 
* @param spaces int     the number of spaces to prepend before printing
*/
void MassDawgNode::show(int spaces){
    for (int i = 0; i < spaces; i++) cout << " ";

    // print kmers
    cout << "|---> kmers: {";
    if (this->kmers.size() > 0) cout << this->kmers[0];
    for (int i = 1; i < this->kmers.size(); i ++) cout << ", " + this->kmers[i];
    cout << "}\n";

    // print the edges
    for (int j = 0; j < this->edges.size(); j ++){
        for (int i = 0; i < spaces + 2; i ++) cout << " ";
        cout << "edge masses: " + to_string(this->edges[j]->singlyMass) + ", " + to_string(this->edges[j]->doublyMass) + "\n";
        // print the child of the edge
        this->edges[j]->child->show(spaces + 4);
    }
    
}
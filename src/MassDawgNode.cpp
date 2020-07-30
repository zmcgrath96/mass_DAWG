#include "MassDawgNode.hpp"
#include "utils.cpp"


MassDawgNode::MassDawgNode (){}

// init with a string
MassDawgNode::MassDawgNode (float singlyMass, float doublyMass, string kmer){
        this->kmers.push_back(kmer);
        this->singlyMass = singlyMass;
        this->doublyMass = doublyMass;
    }

// it is assumed all nodes are deleted INDEPENDENTLY of eachother, 
// nodes are not recursively deleted, so all pointer are set to null
// and that is all
MassDawgNode::~MassDawgNode(){
    try{
        for (MassDawgNode * i: this->children){
        i = nullptr;
    }
    }
    catch (exception) {
        cout << "\nERROR WHEN SETTING POINTERS TO NULL\n";
    }
    
}

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
 * @param singlyMass    float  singly charged mass to connect nodes
 * @param doublyMass    float  doubly charged mass to connect nodes
 * @param kmer          string  the string to associate with the new kmer
 * 
 * @return Edge *       edge connecting the parent to the new child
*/
MassDawgNode * MassDawgNode::addChild(float singlyMass, float doublyMass, string kmer){
    MassDawgNode * newChild = new MassDawgNode(singlyMass, doublyMass, kmer);
    this->children.push_back(newChild);

    return newChild;
}

/**
 * Turn the node's value into a large string in order to make 
 * it hashable
 * 
 * @return string   the hashable string
*/
string MassDawgNode::hash(){
    return to_string(this->singlyMass) + '_' + to_string(this->doublyMass);
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
    cout << "} \t masses: " + to_string(this->singlyMass) + ", " + to_string(this->doublyMass) + "\n";

    // show each child
    for (int i = 0; i < this->children.size(); i++) this->children[i]->show(spaces+2);
}
#include "Edge.hpp"
#include "MassDawgNode.hpp"

// new edge initialized with parameters
Edge::Edge (double singlyMass, double doublyMass, string kmer){
    this->singlyMass = singlyMass;
    this->doublyMass = doublyMass;

    // create a new child node and give it a kmer
    this->child = new MassDawgNode(kmer);
}
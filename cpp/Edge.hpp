#ifndef EDGE_H
#define EDGE_H

#include <string>

class MassDawgNode;

using namespace std;

class Edge {
public:
    double singlyMass;
    double doublyMass;
    MassDawgNode * child;

    // new edge initialized with parameters
    Edge (double singlyMass, double doublyMass, string kmer);
};
#endif
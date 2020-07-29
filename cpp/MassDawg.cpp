#include "MassDawg.hpp"

// empty constructor takes no values
MassDawg::MassDawg(){
    this->root = new MassDawgNode();
}

MassDawg::~MassDawg(){
    delete this->root;
}

/**
 * Show the graph as a tree in the console
*/
void MassDawg::show(){
    cout << "\nroot\n";
    for (int i  = 0; i < this->root->edges.size(); i ++){
        for (int j = 0; j < 2; j++) cout << " ";
        cout << "edge masses: " + to_string(this->root->edges[i]->singlyMass) + ", " + to_string(this->root->edges[i]->doublyMass) + "\n";
        this->root->edges[i]->child->show(4);
    }
}

/**
 * Add a new singly and doubly charged sequence associated with the kmer to the graph
 * 
 * @param singlySequence    vector<doubly>  the singly charged sequence of masses
 * @param doublySequence    vector<doubly>  the doubly charged sequence of masses
 * @param kmer              string          the sequence of amino acids associated with this mass
*/
void MassDawg::insert(vector<double> singlySequence, vector<double> doublySequence, string kmer){
    // if the new one is not greater than the last one, we need to throw
    if (!this->previousIsLessThan(singlySequence, doublySequence)) {
        cout << "Error: sequences must be inserted in lowest to highest order";
        throw;
    }

    // find the common prefix between the new sequence and the last sequence (mass based)
    int commonPrefix = 0;
    int iterLength = MIN(singlySequence.size(), this->ps.singlySequence.size());

    // go through and see how much these sequences have in common
    for (int i = 0; i < iterLength; i++){
        // if either of the singly or doubly sequences are not the same, break
        if ((singlySequence[i] != this->ps.singlySequence[i]) 
        || (doublySequence[i] != this->ps.doublySequence[i])) break;

        // update the kmer at the node at this position in the previous sequence
        (*this->ps.nodes[i]).addKmer(kmer.substr(0, i + 1));

        commonPrefix ++;
    }

    // minimize the previous sequence
    this->minimize(commonPrefix);

    // starting point in the graph for insertion
    MassDawgNode * currentNode;
    // if we don't have any unchecked nodes remaining, set our current 
    // node to root, otherwise point to the last unchecked node
    if (this->uncheckedNodes.size() == 0) currentNode = this->root;
    else currentNode = this->uncheckedNodes.back().child;

    // create a new previous sequence object to keep track of what 
    // we are doing now. We don't need all of the nodes, just the 
    // nodes up until commonPrefix
    PreviousSequence nextPreviousSequence;
    nextPreviousSequence.singlySequence = singlySequence; 
    nextPreviousSequence.doublySequence = doublySequence;

    for (int i = 0; i < commonPrefix; i++){
        nextPreviousSequence.nodes.push_back(this->ps.nodes[i]);
    }

    // go through the remainder of the sequence and create new nodes
    for (int i = commonPrefix; i < singlySequence.size(); i ++){

        // add a new child to my current node
        Edge * newEdge = currentNode->addChild(
            singlySequence[i], 
            doublySequence[i], 
            kmer.substr(0, i + 1)
        );

        // add the pointer to the next node to previous nodes
        nextPreviousSequence.nodes.push_back(newEdge->child);

        // create the another unchecked node and add it to the list
        UncheckedNode un;
        un.child = newEdge->child;
        un.parent = currentNode;
        un.connector = newEdge;
        this->uncheckedNodes.push_back(un);

        currentNode = newEdge->child;
    }

    // update the previous sequence to be this sequence
    this->ps = nextPreviousSequence;
}

/**
 * What makes this a graph and not a tree. Combines nodes that share edges and values
 * 
 * @param downTo    int     the level down to which we need to minimize 
*/
void MassDawg::minimize(int downTo){
    if (this->uncheckedNodes.size() == 0) return;

    // go through all of the unckecked nodes that exist and see if we can merge any
    for (int i = this->uncheckedNodes.size() - 1; i > downTo - 1; i--){

        // local variables to make things easier
        UncheckedNode currentUnchecked = this->uncheckedNodes.back();
        MassDawgNode * child = currentUnchecked.child;
        Edge * parentToChild = currentUnchecked.connector;

        // get the hashable value of the child 
        string childsHash = child->hash();

        // check to see if this value can be found
        unordered_map<string, MassDawgNode *>::const_iterator result = this->minimizedNodes.find(childsHash);
        
        // if the node is not in the map, add it
        if (result == this->minimizedNodes.end()){
            this->minimizedNodes.insert({childsHash, child});
        }

        // Edge * in the current unchecked points to the child. we need to point
        // this to the node in the map (map has a pointer, so just set edge->child to map return)
        // and delete child. Add all the kmers of the child to the one in the map
        else {
            // get the pointer from the map
            MassDawgNode * minimizedNode = result->second;
            // set add all the kmers in the child to the node
            for (int j = 0; j < child->kmers.size(); j++){
                minimizedNode->addKmer(child->kmers[j]);
            }
            // delete child
            delete child;
            // point the edge to the minimized node
            parentToChild->child = minimizedNode;
        }

        // remove the one we just finished working on
        this->uncheckedNodes.pop_back();
    }
}

/**
 * Checks to see if the new sequences are greater than the old previous sequence
 * 
 * @param singlySequence    vector<double>  the new sequence of singly charged masses
 * @param doublySequence    vector<double>  the new sequence of doubly charged masses
 * 
 * @return bool     True if the new sequences are greater than the prvious, False otherwise
*/
bool MassDawg::previousIsLessThan(vector<double> singlySequence, vector<double> doublySequence){
    // get the lengths of each and determine the shorter one
    int newLength = singlySequence.size();
    int oldLength = this->ps.singlySequence.size();

    int iterLength = MIN(newLength, oldLength);

    for (int i = 0; i < iterLength; i ++){

        // if the new one in either the singly OR doubly
        // is greater than the previous at i, retutn False
        if ((singlySequence[i] < this->ps.singlySequence[i])
        || (doublySequence[i] < this->ps.doublySequence[i])){
            return false;
        }
    }

    return newLength >= oldLength;
}


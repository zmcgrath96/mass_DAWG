#include "MassDawgc.hpp"

/**
 * Used for qsort on floats. If d1 < d2, a number < 0 returned
 * if d1 > d2, a number > 0 returned
 * 
*/
int dblCmp(const void * d1, const void * d2){
    return *(float *)d1 - *(float *)d2;
}

/**
 * Convert the ppm value of a given mass to a value in daltons
 * 
 * @param mass      float  the mass to calculate the dalton tolerance for
 * @param ppmTol    int     the tolerance in ppm
 * 
 * @return float   the tolerance in daltons 
*/
float ppmToDa(float mass, int ppmTol){
    return ((float)ppmTol / 1000000.0) * mass;
}


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

/*******************Public methods*******************/

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
    for (int i  = 0; i < (int)this->root->children.size(); i ++) this->root->children[i]->show(2);
}

/**
 * Add a new singly and doubly charged sequence associated with the kmer to the graph
 * 
 * @param singlySequence    vector<doubly>  the singly charged sequence of masses
 * @param doublySequence    vector<doubly>  the doubly charged sequence of masses
 * @param kmer              string          the sequence of amino acids associated with this mass
*/
void MassDawg::insert(vector<float> singlySequence, vector<float> doublySequence, string kmer){
    // if the new one is not greater than the last one, we need to throw
    if (!this->previousIsLessThan(singlySequence, doublySequence)) {
        throw "Error: sequences must be inserted in lowest to highest order";
    }

    // find the common prefix between the new sequence and the last sequence (mass based)
    int commonPrefix = 0;
    int iterLength = MIN((int)singlySequence.size(), (int)this->ps.singlySequence.size());

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
    for (int i = commonPrefix; i < (int)singlySequence.size(); i ++){

        // add a new child to my current node
        MassDawgNode * newChild = currentNode->addChild(
            singlySequence[i], 
            doublySequence[i], 
            kmer.substr(0, i + 1)
        );

        // add the pointer to the next node to previous nodes
        nextPreviousSequence.nodes.push_back(newChild);

        // create the another unchecked node and add it to the list
        UncheckedNode un;
        un.child = newChild;
        un.parent = currentNode;
        this->uncheckedNodes.push_back(un);

        currentNode = newChild;
    }

    // update the previous sequence to be this sequence
    this->ps = nextPreviousSequence;
}

/**
 * Any remaining unchecked nodes will be checked for merging to 
 * complete the dawg. 
*/
void MassDawg::finish(){
    this->minimize(0);
}

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
vector<string> MassDawg::fuzzySearch(vector<float> sequence, int gapAllowance, int ppmTol){
    // save all results from all children into a vector
    vector<vector<string> > allResults;
    for (int i = 0; i < (int)this->root->children.size(); i ++) 
        allResults.push_back(this->fuzzySearchRec(sequence, this->root->children[i], 0, gapAllowance, ppmTol));

    // combine them all into one vector
    vector<string> mergedResults;
    for (int i = 0; i < (int)allResults.size(); i++){
        for (int j = 0; j < (int)allResults[i].size(); j++){
            if (allResults[i][j].empty()) continue;
            mergedResults.push_back(allResults[i][j]);
        } 
    }

    return mergedResults;

}



/*******************Private methods*******************/



/**
 * What makes this a graph and not a tree. Combines nodes that share edges and values
 * 
 * @param downTo    int     the level down to which we need to minimize 
*/
void MassDawg::minimize(int downTo){
    if ((int)this->uncheckedNodes.size() == 0) return;

    // go through all of the unckecked nodes that exist and see if we can merge any
    for (int i = (int)this->uncheckedNodes.size() - 1; i > downTo - 1; i--){

        // local variables to make things easier
        UncheckedNode currentUnchecked = this->uncheckedNodes.back();
        MassDawgNode * child = currentUnchecked.child;
        MassDawgNode * parent = currentUnchecked.parent;

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
            for (int j = 0; j < (int)child->kmers.size(); j++){
                minimizedNode->addKmer(child->kmers[j]);
            }
            
            // find the child of the parent that pointed to the node
            // and update that pointer to point to minimizedNode
            for (int j = 0; j < (int)parent->children.size(); j ++){
                if (parent->children[j] == child) {
                    parent->children[j] = minimizedNode;
                    break;
                }
            }
            
            // delete child
            try {
                delete child;
            }
            catch (exception){
                cout << "Error when deleting child in MassDawg";
            }
        }

        // remove the one we just finished working on
        this->uncheckedNodes.pop_back();
    }
}

/**
 * Checks to see if the new sequences are greater than the old previous sequence
 * 
 * @param singlySequence    vector<float>  the new sequence of singly charged masses
 * @param doublySequence    vector<float>  the new sequence of doubly charged masses
 * 
 * @return bool     True if the new sequences are greater than the prvious, False otherwise
*/
bool MassDawg::previousIsLessThan(vector<float> singlySequence, vector<float> doublySequence){
    // get the lengths of each and determine the shorter one
    int newLength = (int)singlySequence.size();
    int oldLength = (int)this->ps.singlySequence.size();

    int iterLength = MIN(newLength, oldLength);

    for (int i = 0; i < iterLength; i ++){

        // if the new one in either the singly OR doubly
        // is greater than the previous at i, retutn False
        if ((singlySequence[i] < this->ps.singlySequence[i])
        || (doublySequence[i] < this->ps.doublySequence[i])){
            return false;
        }

        // if the new one is greater than the old one, return true
        if ((singlySequence[i] > this->ps.singlySequence[i])
        && (doublySequence[i] > this->ps.doublySequence[i])){
            return true;
        }
    }

    return newLength >= oldLength;
}


/**
 * Recursive search of the graph allowing for gapAllowance missed masses in the
 * search before returning whatever is found at the level
 * 
 * @param sequence      vector<float>  The sequence to use to navigate the graph
 * @param currentNode   MassDawgNode *  The current node to investigate
 * @param currentGap    int             The number of gaps we have allowed up until this point
 * @param gapAllowance  int             The total number of gaps to allow
 * 
 * @return vector<string>   The kmers associated with the deepest part of the branch investigated
*/
vector<string> MassDawg::fuzzySearchRec(vector<float> sequence, MassDawgNode * currentNode, int currentGap, int gapAllowance, int ppmTol){
    // for the cases when we return nothing
    vector<string> emptyResult = {""};

    // BASE CASE: we're past our limit
    if ((gapAllowance - currentGap) < 0) return emptyResult;

    // BASE CASE: we're given an empty sequence
    if (sequence.empty()) return emptyResult;

    // check to see if any of the values in the sequence are within
    // the range of the singly and doubly masses within this node
    float singlyDaTol = ppmToDa(currentNode->singlyMass, ppmTol);
    float doublyDaTol = ppmToDa(currentNode->doublyMass, ppmTol);

    // calcuate the bounds
    float singlyLowerBound = currentNode->singlyMass - singlyDaTol;
    float singlyUpperBound = currentNode->singlyMass + singlyDaTol;
    float doublyLowerBound = currentNode->doublyMass - doublyDaTol;
    float doublyUpperBound = currentNode->doublyMass + doublyDaTol;

    bool massFound = false;
    // go through each mass in the sequence and see if any of the values are in 
    // either set of bounds
    for (int i = 0; i < (int)sequence.size(); i++){
        if ((singlyLowerBound <= sequence[i] && sequence[i] <= singlyUpperBound) 
        || (doublyLowerBound <= sequence[i] && sequence[i] <= doublyUpperBound)){
            massFound = true;
            break;
        }
    }

    // add to the gap if we didnt find the mass
    int gapAddition = massFound ? 0 : 1;

    // updated vector. Won't change if mass wasnt found
    vector<float> updatedSequence;

    // if we found the mass, update sequence to not contain
    // any of the masses < our doubly lower bound and any masses in our singly range
    if (massFound) {
        for (int i = 0; i < (int)sequence.size(); i ++){
            if ((singlyLowerBound <= sequence[i] && sequence[i] <= singlyUpperBound)
             || sequence[i] <= doublyLowerBound){
                 continue;
            }

            //otherwise keep it
            updatedSequence.push_back(sequence[i]);
        }
    }
    else updatedSequence = sequence;

    // if our updated sequence is EMPTY but we found the mass, return my kmers
    if (updatedSequence.empty() and massFound) return vector<string>(currentNode->kmers);

    // otherwise go through all of the children and save their results
    vector<vector<string> > childrensResults;
    for (int i = 0; i < (int)currentNode->children.size(); i++){
        childrensResults.push_back(this->fuzzySearchRec(
            updatedSequence, 
            currentNode->children[i], 
            currentGap + gapAddition, 
            gapAllowance, 
            ppmTol
        ));
    }

    // combine all the children's return values into one vector
    vector<string> results;
    for (int i = 0; i < (int)childrensResults.size(); i++){
        for (int j = 0; j < (int)childrensResults[i].size(); j++){
            if (childrensResults[i][j].empty()) continue;
            results.push_back(childrensResults[i][j]);
        }
    }

    // if we don't have any results and we found a mass, return my results
    if (results.empty() && massFound) return vector<string>(currentNode->kmers);

    // otherwise return results
    return results.empty() ? emptyResult : results;
}

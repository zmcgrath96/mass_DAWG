#include "MassDawg.hpp"
#include "utils.hpp"

/*******************Public methods*******************/

LongestCommonPrefix::LongestCommonPrefix(vector<float> sS, vector<float> dS, vector<MassDawgNode *> nodes){
        this->singlySequence = sS;
        this->doublySequence = dS;
        this->nodes = nodes;
    };

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
    int commonPrefix = 0;
    LongestCommonPrefix lcp;

    // if the new seqeunce is greater than the old sequence, use the this->previousSequence
    // instance with its nodes for speed up in sorted input
    if (this->previousIsLessThan(singlySequence, doublySequence)) {

        // find the common prefix between the new sequence and the last sequence (mass based)	
        int iterLength = MIN((int)singlySequence.size(), (int)this->previousSequence.singlySequence.size());	

        // go through and see how much these sequences have in common	
        for (int i = 0; i < iterLength; i++){	
            // if either of the singly or doubly sequences are not the same, break	
            if ((singlySequence[i] != this->previousSequence.singlySequence[i]) 	
            || (doublySequence[i] != this->previousSequence.doublySequence[i])) break;	

            // update the kmer at the node at this position in the previous sequence
            (*this->previousSequence.nodes[i]).addKmer(kmer.substr(0, i + 1));

            commonPrefix ++;	
        }
        
        // at the break, make lcp the previous sequence
        for (int i = 0; i < commonPrefix; i ++){
            lcp.nodes.push_back(this->previousSequence.nodes[i]);
            lcp.singlySequence.push_back(singlySequence[i]);
            lcp.doublySequence.push_back(doublySequence[i]);
        }
    }	

    // otherise, in the case that the new sequence is smaller than previous (unsorted input)
    // find the longest common prefix
    else {
        // get the longest common prefix of this new sequence
        lcp = this->longestCommonPrefix(singlySequence, doublySequence);

        // add this kmer to all of the nodes in the lcp
        for (MassDawgNode * node: lcp.nodes) node->addKmer(kmer);
        commonPrefix = (int)lcp.singlySequence.size();
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

    // copy over all of the overlapped nodes
    for (MassDawgNode * node: this->previousSequence.nodes){
        nextPreviousSequence.nodes.push_back(node);
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
    this->previousSequence = nextPreviousSequence;
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

/**
* A search with no gaps allowed
* 
* @param sequence       vector<float>   the sequence to search
* @param ppmTol         int             the tolerance in parts per million to accept when searching
* 
* @return vector<string>                All kmers that we found in the search
*/
vector<string> MassDawg::search(vector<float> sequence, int ppmTol){
    MassDawgNode * currentNode = this->root;

    if (sequence.empty()) return vector<string> {};

    // sort the sequence to ensure order
    sort(sequence.begin(), sequence.end(), greater<float>());

    while (true){
        bool extended = false;

        // go through each child and see if the masses fit the tolerance
        for (MassDawgNode * child: currentNode->children){
            // check to see if any of the values in the sequence are within
            // the range of the singly and doubly masses within this node
            float singlyDaTol = ppmToDa(child->singlyMass, ppmTol);
            float doublyDaTol = ppmToDa(child->doublyMass, ppmTol);

            // calcuate the bounds
            float singlyLowerBound = child->singlyMass - singlyDaTol;
            float singlyUpperBound = child->singlyMass + singlyDaTol;
            float doublyLowerBound = child->doublyMass - doublyDaTol;
            float doublyUpperBound = child->doublyMass + doublyDaTol;

            // if any of the masses in the sequence are within this tolerance, we will 
            // continue with this child.
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

            if (!massFound) continue;
            extended = true;

            // reduce our sequence to be those masses that are > doubly upper and not in the singly range
            vector<float> updatedSequence;
            for (float mass: sequence){
                if (mass > doublyUpperBound && 
                !(singlyLowerBound <= mass && mass <= singlyUpperBound)){
                    updatedSequence.push_back(mass);
                }
            }

            sequence = updatedSequence;
            currentNode = child;
            break;
        }

        if (!extended) break;
    }

    return currentNode->kmers;
    
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

        // Parent contains child in one of its children. we need to point
        // this to the node in the map (map has a pointer, so just set parent[child]->child to map return)
        // and delete child. Add all the kmers of the child to the one in the map
        else {
            // get the pointer from the map
            MassDawgNode * minimizedNode = result->second;

            // set add all the kmers in the child to the node
            for (string kmer: child->kmers){
                minimizedNode->addKmer(kmer);
            }
           
            // add all the children of the current child node to the minimized node
            for (MassDawgNode * childsChild: child->children){
                minimizedNode->addChildByPointer(childsChild);
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

    // updated vector. Won't change if mass wasn't found
    vector<float> updatedSequence;

    // if we found the mass, update sequence to not contain
    // any of the masses < our doubly lower bound and any masses in our singly range
    if (massFound) {
        for (int i = 0; i < (int)sequence.size(); i ++){
            if ((singlyLowerBound <= sequence[i] && sequence[i] <= singlyUpperBound) 
            || (doublyLowerBound <= sequence[i] && sequence[i] <= doublyUpperBound)){
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
    int oldLength = (int)this->previousSequence.singlySequence.size();	

    int iterLength = MIN(newLength, oldLength);	

    for (int i = 0; i < iterLength; i ++){	

        // if the new one in either the singly OR doubly	
        // is greater than the previous at i, retutn False	
        if ((singlySequence[i] < this->previousSequence.singlySequence[i])	
        || (doublySequence[i] < this->previousSequence.doublySequence[i])){	
            return false;	
        }	

        // if the new one is greater than the old one, return true	
        if ((singlySequence[i] > this->previousSequence.singlySequence[i])	
        && (doublySequence[i] > this->previousSequence.doublySequence[i])){	
            return true;	
        }	
    }	

    return newLength >= oldLength;	
}

/**
 * Find the longest common prefix of input sequences to a path in the tree. Used for out
 * of order insertions.
 * 
 * @param singlySequence    vector<float>   the singly sequence we are looking for a prefix of
 * @param doublySequence    vector<float>   the doubly sequence we are looking for a prefix of
 * 
 * @returns LongestCommonPrefix *   the class instance holding the longest common prefixe
*/
LongestCommonPrefix MassDawg::longestCommonPrefix(vector<float> singlySequence, vector<float> doublySequence){
    // first check to see if we should return the root
    bool rootHasChild = false;
    float delta = .0001;
    for (MassDawgNode * child: root->children){
        if (abs(child->singlyMass - singlySequence[0]) < delta && abs(child->doublyMass - doublySequence[0]) < delta){
            rootHasChild = true;
            break;
        }
    }

    // if the root does not have a child
    if (!rootHasChild){
        return LongestCommonPrefix(vector<float> {}, vector<float> {}, vector<MassDawgNode *> {this->root});
    }

    MassDawgNode * currentNode = this->root;
    LongestCommonPrefix lcp(vector<float> {}, vector<float> {}, vector<MassDawgNode *> {});
    for (int i = 0; i < (int)singlySequence.size(); i++){
        // go through the current node's children, and if the child has a node that meets the weight at i,
        // add it ot the list

        // keep track if we continued the common prefix
        bool extended = false;

        for (MassDawgNode * child: currentNode->children){

            // if both the singly and doubly sequence were found, add it to the list and end this loop
            if (abs(child->singlyMass - singlySequence[i]) < delta && abs(child->doublyMass - doublySequence[i]) < delta){
                lcp.nodes.push_back(child);
                lcp.singlySequence.push_back(singlySequence[i]);
                lcp.doublySequence.push_back(doublySequence[i]);
                extended = true;
                break;
            }
        }

        // if we didn't extend, break
        if (!extended) break;
    }
    return lcp;
}


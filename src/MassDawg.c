#include <stdio.h>
#include <string.h>

#include "MassDawgNode.c"

// for keeping track of the nodes in the last sequence as well 
// as keeping track of the last sequnece
struct PreviousSequence {
    double * singlySequence;
    double * doublySequence;
    MassDawgNode ** nodes;
    int length; 
};

// for keeping track of nodes that need to be minimized (reduced)
struct UncheckedNode {
    MassDawgNode * parent;
    MassDawgNode * child;
    double singlyEdgeMass; 
};

// for holding all nodes and edges
typedef struct massDawg {
    MassDawgNode * root;
    struct PreviousSequence * previousSequence;
    struct UncheckedNode ** uncheckedNodes;
    int numUncheckedNodes;
    MassDawgNode ** minimizedNodes;
    int numMinimizedNodes;
} MassDawg;

/**
 * Create a new previous sequence struct with the number of 
 * sequences and MassDawgNodes specified
 * 
 * @param length    int     the number of elements of this 
 * 
 * @return PreviousSequence * 
*/
struct PreviousSequence * initPreviousSequence(int length) {
    struct PreviousSequence * ps = malloc(sizeof(struct PreviousSequence));
    MassDawgNode ** nodes = malloc(sizeof(MassDawgNode *) * length);
    double singlySequence[length];
    double doublySequence[length];

    ps->singlySequence = singlySequence;
    ps->doublySequence = doublySequence;
    ps->nodes = nodes;
    ps->length = length;
    
    return ps;
}

/**
 * Allocate space for an unchecked node and set all attributes ot
 * 0s and NULL
 * 
 * @return UncheckedNode *
*/
struct UncheckedNode * initUncheckedNode(){
    struct UncheckedNode * un = malloc(sizeof(struct UncheckedNode));
    un->parent = NULL;
    un->child = NULL;
    un->singlyEdgeMass = 0.0;

    return un;
}

/**
 * Create a new MassDawg
 * 
 * @returns MassDawg *
*/
MassDawg * initMassDawg(){
    MassDawg * md = malloc(sizeof(MassDawg));

    // init the structs 
    md->root = initMassDawgNode();
    md->previousSequence = initPreviousSequence(0);

    // init the pointers to dynamic structs
    md->minimizedNodes = NULL;
    md->uncheckedNodes = NULL;

    // init counters
    md->numUncheckedNodes = 0;
    md->numMinimizedNodes = 0;
    return md;
}

/**
 * Delete a previousSequence struct, pointing all pointers to NULL
 * instead of freeing the nodes. 
 * 
 * @param ps    struct PreviousSequence *   the previous sequence
*/
void clearPreviousSequnce(struct PreviousSequence * ps){
    for (int i = 0; i < ps->length; i ++){
        ps->nodes[i] = NULL;
    }
    free(ps->nodes);
    ps->nodes = NULL;
    free(ps);
    ps = NULL;
}

/**
 * Delete an UncheckedNode struct and set its pointers to NULL
 * 
 * @param un    struct UncheckedNode *  the node to clear
*/
void clearUncheckedNode(struct UncheckedNode * un){
    un->parent = NULL;
    un->child = NULL;
    free(un);
    un = NULL;
}


/**
 * Recursively print the graph as a tree
 * 
 * @param md    MassDawg *  the mass dawg to print
 */
void showDawg(MassDawg * md){
    printf("root\n");
    
    for (int i = 0; i < md->root->numEdges; i ++){
        for (int space = 0; space < 2; space ++) printf(" ");
        printf("edge: {singly: %f, doubly: %f}\n", md->root->edges[i]->singlyMass, md->root->edges[i]->doublyMass);
        showNode(md->root->edges[i]->child, 4);
    }
}

/**
 * Combine all possible nodes that can be combined from here down
 * 
 * @param md        MassDawg *  the graph we are reducing
 * @param downTo    int         the level of the graph down to which we should look to combine nodes
*/
void minimize(MassDawg * md, int downTo){

    for (int i = md->numUncheckedNodes - 1; i > downTo - 1; i--){

        struct UncheckedNode * checkingNode = md->uncheckedNodes[i];

        // if we have minimized the node, set to 1 so that we know to 
        // add this child to the set of minimized 
        int minimized = 0;

        // check to see if this node is in the minimized nodes
        for (int nodeNum = 0; nodeNum < md->numMinimizedNodes; nodeNum ++){
            if (nodesEqual(md->minimizedNodes[nodeNum], checkingNode->child)){
            
                // find the edge with the correct mass we are looking for
                struct Edge * updatingEdge;
                for (int j = 0; j < checkingNode->parent->numEdges; j++){

                    struct Edge * edgeInQuestion = checkingNode->parent->edges[j];

                    if (edgeInQuestion->singlyMass == checkingNode->singlyEdgeMass){
                        // BEFORE we lose the child node, lets merge
                        // add all of the kmers FROM edgeInQuestion's child
                        // TO the minimized node
                        for (int kmerCounter = 0; kmerCounter < edgeInQuestion->child->numKmers; kmerCounter ++){
                            addKmer(md->minimizedNodes[nodeNum], edgeInQuestion->child->kmers[kmerCounter]);
                        }

                        // remove child from the edge
                        deleteMassDawgNode(&edgeInQuestion->child);

                        // update this edge's child node to point to be the minimized node
                        // with this value
                        edgeInQuestion->child = md->minimizedNodes[nodeNum];
                        minimized = 1;

                        break;
                    }
                }

            }

            // if we minimized, break
            if (minimized == 1) break;
        }

        // if we haven't minimized this node yet, add the child to minimized
        if (minimized == 0){

            // if mass minimized nodes haven't been allocated, allocate it
            if (md->minimizedNodes == NULL){
                md->minimizedNodes = malloc(sizeof(MassDawgNode *));
            }
            
            else {
                md->minimizedNodes = realloc(
                    md->minimizedNodes, 
                    sizeof(MassDawgNode *) * (md->numMinimizedNodes + 1)
                );
            }
            
            // append this minmized node to the end of our list
            md->minimizedNodes[md->numMinimizedNodes] = checkingNode->child;
            md->numMinimizedNodes ++;
        }

        // remove the unchecked node from list
        clearUncheckedNode(md->uncheckedNodes[md->numUncheckedNodes-1]);

        // reallocate memory such that the we have 1 less node
        md->uncheckedNodes = realloc(
            md->uncheckedNodes, 
            sizeof(struct UncheckedNode *) * (md->numUncheckedNodes - 1)
        );
    
        md->numUncheckedNodes --;
    }
}

/**
 * Insert a sequence of masses into the Dawg. The input 
 * to the insertion (for many insertions) must be sorted.
 * 
 * @param md                MassDawg * the MassDawg struct to insert into
 * @param singlySeqeunce    double * the list of doubles that make up the singly charged spectrum
 * @param doublySequence    double * the list of doubles that make up the doubly charged spectrum
 * @param kmer              char * the sequence of amino acids associated with these masses
 * @param sequenceLength    int the length of the sequence
*/
void insert(MassDawg * md, double * singlySequence, double * doublySequence, char * kmer, int sequenceLength){

    // if the previous sequence is not none (from init) and the 
    // new sequence is > the old sequence, then we need to error
    if (md->previousSequence->length > 0 
        && (compareDoubleArrays(
                singlySequence, 
                md->previousSequence->singlySequence, 
                sequenceLength, 
                md->previousSequence->length
            )
        || compareDoubleArrays(
                doublySequence, 
                md->previousSequence->doublySequence,
                sequenceLength, 
                md->previousSequence->length
    ))){
        printf("ERROR: Input sequence must be less than the last sequence.");
        return;
    }

    // Continue with adding the sequnce since the sequence is greater than the last

    // keep track of how long the common prefix is between the old sequence and
    // the new sequence
    int commonPrefixLength = 0;

    // minimum iterative length between the two sequences
    int iterLength = MIN(
        sequenceLength, 
        md->previousSequence->length
    );

    // go through each of the values in the two sequences and compare values until
    // we either run out of values or the two differ
    for (int i = 0; i < iterLength; i++) {
        
        // check to see if values at the two sequences at position i are the same.
        // if not, break 
        if (singlySequence[i] != md->previousSequence->singlySequence[i]
            || doublySequence[i] != md->previousSequence->doublySequence[i]) break;

        // add the current kmer from 0 to i to the node at i in the previous sequence
        // start by getting the correct kmer from the string passed in
        char* newKmer = malloc((i+2) * sizeof(char));
        strncpy(newKmer, kmer, i+1);
        newKmer[i+1] = '\0'; 

        // add the kmer to the node. MassDawgNode checks for duplicates, so we just need to 
        // pass the kmer in here and let it handle the rest
        addKmer(md->previousSequence->nodes[i], newKmer);

        // we've made it this far so we increment
        commonPrefixLength ++;
    }

    // combine all nodes we can combine
    minimize(md, commonPrefixLength);

    // the current node we need to start at so that we can add the suffix without repeats
    MassDawgNode * currentNode;

    //add the suffix, starting from the correct node mid-way through the graph
    if (md->numUncheckedNodes == 0) currentNode = md->root;
    // get the child of the last unchecked nodes
    else currentNode = md->uncheckedNodes[md->numUncheckedNodes - 1]->child;

    // create the new previous sequence that will contain all the new information for this new 
    // singly and doubly sequence
    struct PreviousSequence * nextPreviousSequence = initPreviousSequence(sequenceLength);
    nextPreviousSequence->singlySequence = singlySequence;
    nextPreviousSequence->doublySequence = doublySequence;

    // copy over all of the old nodes from the old previous to the new
    for (int i = 0; i < commonPrefixLength; i++){
        nextPreviousSequence->nodes[i] = md->previousSequence->nodes[i];
    }

    clearPreviousSequnce(md->previousSequence);

    for (int i = commonPrefixLength; i < sequenceLength; i++){

        // get the kmer for this set of masses
        char* newKmer = malloc((i+2) * sizeof(char));
        strncpy(newKmer, kmer, i+1);
        newKmer[i+1] = '\0'; 

        // create a new edge and node onto the current node
        MassDawgNode * newChild = addChild(currentNode, newKmer, singlySequence[i], doublySequence[i]);

        // append it to the previous sequence
        nextPreviousSequence->nodes[i] = newChild;

        // create a new unchecked node struct to append to the unckecked nodes 
        struct UncheckedNode * thisUncheckedNode = initUncheckedNode();
        thisUncheckedNode->child = newChild;
        thisUncheckedNode->parent = currentNode;
        thisUncheckedNode->singlyEdgeMass = singlySequence[i];

        // increment the size of the unchecked node list
        md->uncheckedNodes = realloc(
            md->uncheckedNodes, 
            sizeof(struct UncheckedNode *) * (md->numUncheckedNodes + 1)
        );

        // append the new unchecked node to the end of the list and increment our count
        md->uncheckedNodes[md->numUncheckedNodes] = thisUncheckedNode;
        md->numUncheckedNodes ++;

        // set node to the child node
        currentNode = newChild;
    }

    // set the last node's final value to true and the previous sequence to the next previous sequence
    currentNode->final = 1;
    md->previousSequence = nextPreviousSequence;
}

/**
 * Recursively search all child nodes for the input sequence
 * allowing for up to gapLimit missed matches
 * 
 * @param mdn           MassDawgNode *  the current node in our search through the graph
 * @param sequence      double *        the double array we are searching
 * @param sequenceLen
 * @param currentGap    int             the current number of missed sequences in our search
 * @param gapLimit      int             the limit on the number of allowed misses in the search
 * @param ppmTol        int             the tolerance to allow when searching for matches
 * 
 * @return char **                      the list of kmers that exist at the node at the end of the search
*/
char ** __fuzzySearchRec(MassDawgNode * mdn, double * sequence, int sequenceLength, int currentGap, int gapLimit, int ppmTol){
    // place holder for cases when we need to return nothing
    char ** emptyList = NULL;

    // if we've hit our gap limit, return an empty list
    if (gapLimit - currentGap < 0) return emptyList;

    // keep track of edges that fit our criteria
    struct Edge ** edgesWithMass = NULL;
    int * edgesWithMassIndices = NULL;
    int numEdges = 0;

    // the return value
    char ** returnKmers = NULL;
    
    // get all edges that have a mass within our tolerance 
    for (int i = 0; i < sequenceLength; i ++){
        double mass = sequence[i];
        double daTol = ppmToDal(mass, ppmTol);
        double lowerBound = mass - daTol;
        double upperBound = mass + daTol;

        for (int j = 0; j < mdn->numEdges; j++){
            // see if either of the masses (either singly or doubly) are within our bounds
            if ((mdn->edges[j]->singlyMass >= lowerBound && mdn->edges[j]->singlyMass <= upperBound) 
            || (mdn->edges[j]->doublyMass >= lowerBound && mdn->edges[j]->doublyMass <= upperBound)){
                // if we haven't allocated memory for edges with mass, do so
                if (edgesWithMass == NULL) {
                    edgesWithMass = malloc(sizeof(struct Edge *));
                    edgesWithMassIndices = malloc(sizeof(int));
                }
                else {
                    edgesWithMass = realloc(edgesWithMass, sizeof(struct Edge *) * (numEdges + 1));
                    edgesWithMassIndices = realloc(edgesWithMassIndices, sizeof(int) * (numEdges + 1));
                }

                edgesWithMass[numEdges] = mdn->edges[j];
                edgesWithMassIndices[numEdges] = j;
                numEdges ++;
            }
        }
        
    }

    // if we found any edges that met the criteria, go through those and call 
    // fuzzySearchRec again with no increase in the current gap
    if (numEdges > 0){

        // iterate through each edge
        for (int i = 0; i < numEdges; i++){

            // create a new sequence double * that has elements > than the singly mass
            // of the edge and does not contain the doubly mass
            double * newSequence = NULL;
            int newSequenceLength = 0;
            for (int j = 0; j < sequenceLength; j ++){
                // values for checking to see if we should skip this mass
                double edgeSinglyMassUpperBound = edgesWithMass[i]->singlyMass + ppmToDal(edgesWithMass[i]->singlyMass, ppmTol);
                double edgeDoublyMassUpperBound = edgesWithMass[i]->doublyMass + ppmToDal(edgesWithMass[i]->doublyMass, ppmTol);
                double edgeDoublyMassLowerBound = edgesWithMass[i]->doublyMass - ppmToDal(edgesWithMass[i]->doublyMass, ppmTol);
                
                // if the value of sequence at j is <= the singly mass upper bound OR
                // its in the range of the doubly mass, don't add it
                if (edgeSinglyMassUpperBound >= sequence[j]
                || (sequence[j] >= edgeDoublyMassLowerBound && sequence[j] <= edgeDoublyMassUpperBound)){
                    continue;
                }

                // if new sequence is still null, allocate it, otherwise reallocate it
                if (newSequence == NULL) newSequence = malloc(sizeof(double));
                else newSequence = realloc(newSequence, sizeof(double) * (newSequenceLength + 1));
                
                newSequence[newSequenceLength] = sequence[j];
                newSequenceLength ++;
            }

            // get the new return value of the next level recursion
            char ** nextRecursionLevelReturn = __fuzzySearchRec(
                edgesWithMass[i]->child, 
                newSequence, 
                newSequenceLength,
                currentGap, 
                gapLimit, 
                ppmTol
            );

            // if return kmers is NULL, point it to the next level recursion and continue
            returnKmers = concatStrArrays(returnKmers, nextRecursionLevelReturn);
        }
    }

    // for those edges that we did not find a mass, recurse another level with incremented gap
    for (int i = 0; i < mdn->numEdges; i++){
        // check to see if i is in the list of edges that were found to have the correct mass
        int iInEdgeFound = 0;
        for (int j = 0; j < numEdges; j++){
            if (i == edgesWithMassIndices[j]){
                iInEdgeFound = 1;
                break;
            }
        }
        if (iInEdgeFound == 1) continue;

        // get the next level recursion return value
        char ** nextRecursionLevelReturn = __fuzzySearchRec(
            mdn->edges[i]->child, 
            sequence, 
            sequenceLength,
            currentGap + 1, 
            gapLimit, 
            ppmTol
        );

        // concat the two lists
        returnKmers = concatStrArrays(returnKmers, nextRecursionLevelReturn);
    }

    //if return masses is empty but at this depth we found something, return the value of the
    //child of the edges that was found
    if (returnKmers == NULL){
        // go through all edges that had something and append to list
        for (int i = 0; i < mdn->numEdges; i ++){
            char ** kmersCopy = deepCopyStrArray(mdn->edges[i]->child->kmers);
            returnKmers = concatStrArrays(returnKmers, kmersCopy);
        }
    }

    return returnKmers;
}


/**
 * Use a sequence as a guide through the graph. Up to gap
 * missing links are allowed in the graph before returning the
 * value found at that depth. 
 * 
 * @example
 *      input sequence: [100, 200, 500]
 *      graph:          [(100, A), (200, AB), (300, ABC ABX), 
 *                      (400, ABCD ABXY), (500, ABCDE ABXYZ)]
 *      if gap is 0 or 1, the return will be AB
 *      if gap is 2 +, output will be ABCDE ABXYZ
 * 
 * @param md        MassDawg *  the graph to search
 * @param sequence  double *    list of doubles to use to traverse
 * @param gap       int         the number of missing links to allow
 * @param ppmTol    int         the tolerance to allow when searching for matches
 * 
 * @return char **              list of kmers that match the input sequence
*/
char ** fuzzySearch(MassDawg * md, double * sequence, int sequenceLength, int gap, int ppmTol){
    return __fuzzySearchRec(md->root, sequence, sequenceLength, 0, gap, ppmTol);
}

/**
 * Complete the DAWG by minimizing all unchecked nodes that remain. 
 * All minimized and unchecked nodes are cleared
 * 
 * @param md    MassDawg *  the dawg to finish
*/
void finish(MassDawg * md){
    minimize(md, 0);
    
    // clear the unchecked nodes 
    for (int i = 0; i < md->numUncheckedNodes; i++){
        clearUncheckedNode(md->uncheckedNodes[i]);
    }
    free(md->uncheckedNodes);
    md->uncheckedNodes = NULL;
    md->numUncheckedNodes = 0;

}

/**
 * Delete the nodes from the tree and free up all memory
 * 
 * @param md    MassDawg *  the mass dawg to delete
*/
void deleteMassDawg(MassDawg * md){
 
    // all nodes exist in minimized nodes after finish(), so we should
    // call finish on md to ensure that all nodes are in this set
   finish(md);

    //free up all minimized nodes
    for (int i = 0; i < md->numMinimizedNodes; i++){

        // point child of the edges of this node to null so that
        // deleting an edge doesn't recuresively try and delete more
        for (int j = 0; j < md->minimizedNodes[i]->numEdges; j++){
            md->minimizedNodes[i]->edges[j]->child = NULL;
            deleteEdge(md->minimizedNodes[i]->edges[j]);
        }
        deleteMassDawgNode(&md->minimizedNodes[i]);
    }
    free(md->minimizedNodes);
    md->minimizedNodes = NULL;

    // free the previous sequence
    for (int i = 0; i < md->previousSequence->length; i++){
        // all nodes should be deleted, so we just need to 
        // set pointers to null
        md->previousSequence->nodes[i] = NULL;
    }
    free(md->previousSequence->nodes);
    md->previousSequence = NULL;

    // set the root edge pointers to null (the child of those edges must 
    // have been in one of the lists) and delete the root
    for (int i = 0; i < md->root->numEdges; i++){
        md->root->edges[i]->child = NULL;
        deleteEdge(md->root->edges[i]);
    }
    free(md->root->edges);
    md->root->edges = NULL;
    deleteMassDawgNode(&md->root);

    // finally free the struct itself
    free(md);
    md = NULL;
}
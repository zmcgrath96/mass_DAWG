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
struct uncheckedNode {
    MassDawgNode * parent;
    MassDawgNode * child;
    double singlyEdgeMass; 
};

// for holding all nodes and edges
typedef struct massDawg {
    MassDawgNode * root;
    struct PreviousSequence * previousSequence;
    struct uncheckedNode ** uncheckedNodes;
    int numUncheckedNodes;
    MassDawgNode ** minimizedNodes;
    int numMinimizedNodes;
} MassDawg;

/**
 * Create a new previous sequence struct with the number of 
 * sequences and MassDawgNodes specified
 * 
 * @param length int the number of elements of this 
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
 * Delete a previousSequence struct, pointing all pointers to NULL
 * instead of freeing the nodes. 
 * 
 * @param ps struct PreviousSequence * the previous sequence
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
 * Recursively print the graph as a tree
 * 
 * @param md MassDawg * the mass dawg to print
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
 * Recursively delete the nodes from the tree and free up all memory
 * 
 * @param md MassDawg * the mass dawg to delete
*/
void deleteMassDawg(MassDawg * md){
    // first free up all nodes
    deleteMassDawgNode(md->root);

    // free up all unchecked nodes
    for (int i = 0; i < md->numUncheckedNodes; i++){

        // make sure all child and parent nodes are freed
        if (md->uncheckedNodes[i] != NULL){
            if (md->uncheckedNodes[i]->parent != NULL) deleteMassDawgNode(md->uncheckedNodes[i]->parent);
            if (md->uncheckedNodes[i]->child != NULL) deleteMassDawgNode(md->uncheckedNodes[i]->child);
        }
        free(md->uncheckedNodes[i]);
        md->uncheckedNodes = NULL;
    }
    free(md->uncheckedNodes);
    md->uncheckedNodes = NULL;

    //free up all minimized nodes
    for (int i = 0; i < md->numMinimizedNodes; i++){

        // make sure all these nodes are free
        if (md->minimizedNodes[i] != NULL) deleteMassDawgNode(md->minimizedNodes[i]);
    }
    free(md->minimizedNodes);
    md->minimizedNodes = NULL;

    // free the previous sequence
    for (int i = 0; i < md->previousSequence->length; i++){
        if (md->previousSequence->nodes[i] != NULL) deleteMassDawgNode(md->previousSequence->nodes[i]);
    }
    free(md->previousSequence->nodes);
    md->previousSequence = NULL;

    // finally free the struct itself
    free(md);
    md = NULL;
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
 * Combine all possible nodes that can be combined from here down
 * 
 * @param md MassDawg * the graph we are reducing
 * @param downTo int the level of the graph down to which we should look to combine nodes
*/
void minimize(MassDawg * md, int downTo){

    for (int i = md->numUncheckedNodes - 1; i > downTo - 1; i--){

        struct uncheckedNode * checkingNode = md->uncheckedNodes[i];

        // if we have minimized the node, set to 1 so that we know to 
        // add this child to the set of minimized 
        int minimized = 0;

        // check to see if this node is in the minimized nodes
        for (int nodeNum = 0; nodeNum < md->numMinimizedNodes; nodeNum ++){
            if (nodesEqual(md->minimizedNodes[nodeNum], checkingNode->parent)){
            
                // find the edge with the correct mass we are looking for
                struct Edge * updatingEdge;
                for (int j = 0; j < checkingNode->parent->numEdges; j++){
                    struct Edge * edgeInQuestion = checkingNode ->parent->edges[j];
                    if (edgeInQuestion->singlyMass == checkingNode->singlyEdgeMass){
                        // update this edge's child node to be the minimized node
                        // with this value
                        edgeInQuestion->child = md->minimizedNodes[nodeNum];
                        minimized = 1;
                        // remove the child node from unchecked nodes
                        deleteMassDawgNode(checkingNode->child);
                        break;
                    }
                }

            }

            // if we minimized, break
            if (minimized == 1) break;
        }

        // if we haven't minimized this node yet, add the child to minimized
        if (minimized == 0){
            md->numMinimizedNodes ++;

            MassDawgNode ** updatedMinimizedNodes;

            // if mass minimized nodes haven't been allocated, allocate it
            if (md->minimizedNodes == NULL){
                updatedMinimizedNodes = malloc(
                   sizeof(MassDawgNode *)
                );
            }
            
            else {
                printf("Reallocating for updatedMinimized nodes\n");
                updatedMinimizedNodes = realloc(
                    md->minimizedNodes, 
                    sizeof(md->minimizedNodes) + sizeof(MassDawgNode *)
                );
            }
            
            updatedMinimizedNodes[md->numMinimizedNodes-1] = checkingNode->child;
            md->minimizedNodes = updatedMinimizedNodes;
        }

        // remove the unchecked node from list
        free(md->uncheckedNodes[md->numUncheckedNodes - 1]);

        printf("Reallocing for updated unchecked nodes..\n");
        struct uncheckedNode ** updatedUncheckedNodes= realloc(
            md->uncheckedNodes, 
            sizeof(md->uncheckedNodes) - sizeof(struct uncheckedNode *)
        );
        md->numUncheckedNodes --;
        md->uncheckedNodes = updatedUncheckedNodes;
    }
}

/**
 * Insert a sequence of masses into the Dawg. The input 
 * to the insertion (for many insertions) must be sorted.
 * 
 * @param md MassDawg * the MassDawg struct to insert into
 * @param singlySeqeunce double * the list of doubles that make up the singly charged spectrum
 * @param doublySequence double * the list of doubles that make up the doubly charged spectrum
 * @param kmer char * the sequence of amino acids associated with these masses
 * @param sequenceLength int the length of the sequence
*/
void insert(MassDawg * md, double * singlySequence, double * doublySequence, char * kmer, int sequenceLength){

    // if the previous sequence is not none (from init) and the 
    // new sequence is > the old sequence, then we need to error
    if (md->previousSequence->length > 0 
        &&(compareDoubleArrays(
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
        struct uncheckedNode * thisUncheckedNode = malloc(sizeof(struct uncheckedNode));
        thisUncheckedNode->child = newChild;
        thisUncheckedNode->parent = currentNode;
        thisUncheckedNode->singlyEdgeMass = singlySequence[i];

        // increment the size of the unchecked node list
        md->uncheckedNodes = realloc(
            md->uncheckedNodes, 
            sizeof(struct uncheckedNode *) * (md->numUncheckedNodes + 1)
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

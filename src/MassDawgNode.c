#include <string.h>
#include <stdlib.h>

#include "utils.c"

typedef struct massDawgNode {
    struct Edge ** edges;
    int numEdges;
    char ** kmers;
    int numKmers;
    int final;
} MassDawgNode;

struct Edge {
    MassDawgNode * child;
    double singlyMass;
    double doublyMass;
};

/**
 * Allocate memory for a new MassDawgNode and init to 0s and NULL
 * 
 * @return  MassDawgNode*   pointer to the new mass node
 * 
*/
MassDawgNode * initMassDawgNode() {
    MassDawgNode * mdn = malloc(sizeof(MassDawgNode));
    mdn->final = 0;
    mdn->numEdges = 0;
    mdn->numKmers = 0;

    mdn->edges = NULL;
    mdn->kmers = NULL;
    return mdn;
}

/**
 * Allocate memory for a new Edge and init to 0s and NULL
 * 
 * @return  struct Edge *
*/
struct Edge * initEdge(){
    struct Edge * edge = malloc(sizeof(struct Edge));
    edge->doublyMass = 0.0;
    edge->singlyMass = 0.0;
    edge->child = NULL;

    return edge;
}

/**
 * Recursively delete a MassDawgNode and all of its children
 * 
 * @param mdn   MassDawgNode *  the mass dawg node to free
 */
void deleteMassDawgNode(MassDawgNode ** mdn) {
    // delete any and all children
    if ((*mdn)->numEdges > 0){
        // go through all of the edges and call delete on their children        
        for (int i = 0; i < (*mdn)->numEdges; i++){

            // check first to see if someone else deleted this child
            if ((*mdn)->edges == NULL || &(*mdn)->edges[i] == NULL 
                || (*mdn)->edges[i]->child==NULL) continue;
            deleteMassDawgNode(&(*mdn)->edges[i]->child);
            free((*mdn)->edges[i]);
            (*mdn)->edges[i] = NULL;
        }
    }
    if ((*mdn)->edges != NULL){
        free((*mdn)->edges);
        (*mdn)->edges = NULL;
    }
    

    // delete all kmers
    if ((*mdn)->numKmers > 0){
        for (int i = 0; i < (*mdn)->numKmers; i++){
            free((*mdn)->kmers[i]);
            (*mdn)->kmers[i] = NULL;
        }
    }
    if ((*mdn)->kmers != NULL){
        free((*mdn)->kmers);
        (*mdn)->kmers = NULL;
    }
    

    // set ints to 0 just in case
    (*mdn)->numEdges = 0;
    (*mdn)->numKmers = 0;
    (*mdn)->final = 0;

    if(*mdn != NULL) free(*mdn);
    *mdn = NULL;
}

/**
 * Recuresively delete an edge and its children
 * 
 * @param edge  struct Edge *   the edge to delete
*/
void deleteEdge(struct Edge * edge){
    // delete its child first
    if (edge->child != NULL){
        deleteMassDawgNode(&edge->child);;
        edge->child = NULL;
    }

    free(edge);
    edge = NULL;    
}

/**
 * Remove the edge at the provided index from the node
 * 
 * @param mdn       MassDawgNode *  the node to remove the edge from
 * @param edgeIndex int             the index in the edge array to remove
*/
void removeEdgeFromNode(MassDawgNode * mdn, int edgeIndex){
    deleteEdge(mdn->edges[edgeIndex]);
    mdn->edges[edgeIndex] = NULL;

    // shift all of the elements over to the left one index
    for (int shiftedIndex = edgeIndex; shiftedIndex < mdn->numEdges - 1; shiftedIndex ++){
        mdn->edges[shiftedIndex] = mdn->edges[shiftedIndex + 1];
    }
    
    // resize the edges array to be 1 less
    mdn->edges = realloc(
        mdn->edges, 
        sizeof(struct Edge *) * (mdn->numEdges - 1)
    );

    // decrement the edge counter
    mdn->numEdges --;
}

/**
 * Add a new kmer to the list of kmers in a mass dawg node. 
 * Node is checked for the value being passed in to ensure that
 * a kmer is not duplicated
 * 
 * @param mdn   MassDawgNode *  the node to add a kmer to
 * @param kmer  char *          the new string to add 
 * 
 * @return void
 * 
*/
void addKmer(MassDawgNode * mdn, char * kmer){
    // First check to see if this kmer exists in the nodes kmer set
    int kmerFound = 0;
    for (int i = 0; i < mdn->numKmers; i++){
        // check to see if the string is found. 0 is returned from strcmp if they are equal
        if (strcmp(kmer, mdn->kmers[i]) == 0){
            kmerFound = 1;
            break;
        }
    }

    // if the kmer was found, return
    if (kmerFound == 1) return;

    // otherwise create new room for this kmer and add it to the list
    //reallocate the memory for kmers to add a new one
    mdn->kmers = realloc(mdn->kmers, sizeof(char *) * (mdn->numKmers +1));
    char * addingKmer = malloc(sizeof(char) * (strlen(kmer) + 1));
    strcpy(addingKmer, kmer);
    addingKmer[strlen(kmer)] = '\0';
    mdn->kmers[mdn->numKmers] = addingKmer;
    mdn->numKmers ++;
}

/**
 * Add a new edge with a kmer and mass to the node
 * 
 * @param mdn   MassDawgNode *  the node to add an edge to
 * @param kmer  char *          the new kmer associated with the new node
 * @param mass  double          the new mass to associate with the new node
 * 
 * @return MassDawgNode *           the new child node of the new edge
*/
MassDawgNode * addChild(MassDawgNode * mdn, char * kmer, double singlyMass, double doublyMass){

    // create the new mass node
    MassDawgNode * newNode = initMassDawgNode();
    newNode->kmers = malloc(sizeof(char*));
    newNode->kmers[0] = kmer;
    newNode->numKmers ++;

    // make a new edge
    struct Edge * newEdge = initEdge();
    newEdge->singlyMass = singlyMass;
    newEdge->doublyMass = doublyMass;
    newEdge->child = newNode;

    // add the edge to the current node
    // if we don't have any edges yet, malloc memory
    if (mdn->numEdges == 0){
        mdn->edges = malloc(sizeof(struct Edge *));
    }
    else {
        mdn->edges = realloc(mdn->edges, sizeof(mdn->edges) + sizeof(struct Edge *));
    }
    
    mdn->edges[mdn->numEdges] = newEdge;
    mdn->numEdges ++;

    return newNode;
}

/**
 * Determine if a node has a child with a value equal to the input value
 * 
 * @param mdn   MassDawgNode *  the node to check children of 
 * @param mass  double          the mass to check for
 * 
 * @return int 0 if the mass is not found, 1 otherwise
 */ 
int hasChildWithValue(MassDawgNode * mdn, double mass) {
    for (int i = 0; i < mdn->numEdges; i ++){
        if (mass == mdn->edges[i]->singlyMass || mass == mdn->edges[i]->doublyMass){
            return 1;
        }
    }
    return 0;
}

/**
 * Return a child with a mass that child exists. Returns MassDawgNode *. 
 * An empty MassDawgNode pointer is returned if 
 * the child does not exist
 * 
 * @param mdn   MassDawgNode *  the node to check the children of
 * @param mass  double          the mass to retrieve the child of
 * 
 * @return MassDawgNode *
*/
MassDawgNode * getChildWithValue(MassDawgNode * mdn, double mass){
    for (int i = 0; i < mdn->numEdges; i ++){
        if (mass == mdn->edges[i]->singlyMass || mass == mdn->edges[i]->doublyMass){
            return mdn->edges[i]->child;
        }
    }

    return initMassDawgNode();
}

/**
 * Determine if two edges are equl by their values
 * 
 * @param edge1 Edge *  the first edge in the comparison
 * @param edge2 Edge *  the second edge in the comparison
 * 
 * @return int 0 if the two are not equivalent, 1 otherwize
*/
int edgesEqual(struct Edge * edge1, struct Edge * edge2){
    if (edge1->singlyMass == edge2->singlyMass && edge1->doublyMass == edge2->doublyMass){
        return 1;
    }
    return 0;
}

/**
 * Determine if two nodes are equal by their edge values, not kmer values
 * 
 * @param node1 MassDawgNode *  the first node in the comparison
 * @param node2 MassDawgNode *  the second node in the comparison
 * 
 * @return int 0 if the two are not equivalent, 1 otherwise
*/
int nodesEqual(MassDawgNode * node1, MassDawgNode * node2){

    // check the final values first
    if (node1->final != node2->final) return 0;

    // check to see if the num of edges are the same and that each 
    // edge has the same value
    if (node1->numEdges != node2->numEdges) return 0;

    // iterate through all edges to see if all edges are there
    for (int i = 0; i < node1->numEdges; i ++){

        int edgeFound = 0;
        // see if any of the edges are the same
        for (int j = 0; j < node2->numEdges; j ++){
            if (edgesEqual(node1->edges[i], node2->edges[j]) == 1){
                edgeFound = 1;
                break;
            }
        }

        // if the edge wasn't found, return 0
        if (edgeFound == 0) return 0;
    }

    return 1;
}

/**
 * Recursively show this node and all of its descendants
 * 
 * @param mdn       MassDawgNode *  the node to show and its descendants
 * @param spaces    int             the number of spaces to print before
*/
void showNode(MassDawgNode * mdn, int spaces){
    for (int i = 0; i < spaces; i ++) printf(" ");
    printf("|---> k-mers: [");
    
    // print the first one for pretty printing 
    if (mdn->numKmers > 0) printf("%s", mdn->kmers[0]);
    for (int kmerCount = 1; kmerCount < mdn->numKmers; kmerCount ++){
        printf(", %s", mdn->kmers[kmerCount]);
    }
    printf("]\n");

    for (int edgeCount = 0; edgeCount < mdn->numEdges; edgeCount ++){
        for (int i = 0; i < spaces + 2; i ++) printf(" ");

        printf("edge: {singly: %f, doubly: %f}\n", 
            mdn->edges[edgeCount]->singlyMass, 
            mdn->edges[edgeCount]->doublyMass
        );
        showNode(mdn->edges[edgeCount]->child, spaces + 4);
        
    }

}
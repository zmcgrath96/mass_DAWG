#include <string.h>

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

#define EDGE_ARR_SIZE(edges)    ((int)sizeof(edges)/sizeof(struct Edge *))

/**
 * Allocate memory for a new MassDawgNode and init final to 0
 * 
 * @return MassDawgNode* to the new mass node
 * 
*/
MassDawgNode * initMassDawgNode() {
    MassDawgNode * mdn;
    mdn = malloc(sizeof(MassDawgNode));
    mdn->final = 0;
    mdn->numEdges = 0;
    mdn->numKmers = 0;
    return mdn;
}

/**
 * Recursively delete a MassDawgNode and all of its children
 * 
 * @param mdn (MassDawgNode *) the mass dawg node to free
 * 
 * @return void
 * 
 */
void deleteMassDawgNode(MassDawgNode * mdn) {
    // delete any and all children
    if (mdn->numEdges > 0){
        // go through all of the edges and call delete on their children        
        for (int i = 0; i < mdn->numEdges; i++){
            deleteMassDawgNode(mdn->edges[i]->child);
            free(mdn->edges[i]);
            mdn->edges[i] = NULL;
        }
    }
    mdn->edges = NULL;

    // delete all kmers
    if (mdn->numKmers > 0){
        for (int i = 0; i < mdn->numKmers; i++){
            free(mdn->kmers[i]);
            mdn->kmers[i] = NULL;
        }
    }
    mdn->kmers = NULL;

    free(mdn);
    mdn = NULL;
}

/**
 * Add a new kmer to the list of kmers in a mass dawg node
 * 
 * @param kmer char *      the new string to add 
 * @param mdn MassDawgNode *   the node to add a kmer too
 * 
 * @return void
 * 
*/
void addKmer(char * kmer, MassDawgNode * mdn){
    //reallocate the memory for kmers to add a new one
    char ** newKmers = realloc(mdn->kmers, sizeof(mdn->kmers) + sizeof(char*));
    
    // add the new kmers to the end of the kmers
    int numKmers = sizeof(newKmers)/sizeof(char *);
    newKmers[numKmers-1] = kmer;

    // add it to the node
    mdn->kmers = newKmers;
    mdn->numKmers ++;
}

/**
 * Add a new edge with a kmer and mass to the node
 * 
 * @param mdn MassDawgNode *         the node to add an edge to
 * @param kmer char *               the new kmer associated with the new node
 * @param mass double                the new mass to associate with the new node
 * 
 * @return MassDawgNode *           the new child node of the new edge
*/
MassDawgNode * addNewNode(MassDawgNode * mdn, char * kmer, double singlyMass, double doublyMass){

    // create the new mass node
    MassDawgNode * newNode = initMassDawgNode();
    newNode->kmers = malloc(sizeof(char*));
    newNode->kmers[0] = kmer;
    newNode->numKmers ++;

    // make a new edge
    struct Edge * newEdge = malloc(sizeof(struct Edge));
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
 * @param mdn MassDawgNode * the node to check children of 
 * @param mass double the mass to check for
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
 * @param mdn MassDawgNode * the node to check the children of
 * @param mass the mass to retrieve the child of
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
 * @param edge1 Edge * the first edge in the comparison
 * @param edge2 Edge * the second edge in the comparison
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
 * Determine if two nodes are equal by thier values
 * 
 * @param node1 MassDawgNode * the first node in the comparison
 * @param node2 MassDawgNode * the second node in the comparison
 * 
 * @return int 0 if the two are not equivalent, 1 otherwise
*/
int nodesEqual(MassDawgNode * node1, MassDawgNode * node2){

    // check the final values first
    if (node1->final != node2->final) return 0;

    // check each kmer
    int numNode1Kmers = STR_ARR_LEN(node1->kmers);
    int numNode2Kmers = STR_ARR_LEN(node2->kmers);

    if (numNode1Kmers != numNode2Kmers)return 0;

    // check to see if each kmer value is equal
    for (int i = 0; i < numNode2Kmers; i++){

        int strFound = 0;
        // look through all of node 1 kmers to see if they are the same
        for (int j = 0; j < numNode1Kmers; j++){
            if (strcmp(node1->kmers[j], node2->kmers[i]) == 0){
                strFound = 1;
                break;
            }
        }
        
        // if not found, return 0
        if (strFound == 0) return 0;
    }

    // check to see if the num of edges are the same and that each 
    // edge has the same value
    int numNode1Edges = EDGE_ARR_SIZE(node1->edges);
    int numNode2Edges = EDGE_ARR_SIZE(node2->edges);

    if (numNode1Edges != numNode2Edges) return 0;

    // iterate through all edges to see if all edges are there
    for (int i = 0; i < numNode2Edges; i ++){

        int edgeFound = 0;
        // see if any of the edges are the same
        for (int j = 0; j < numNode2Edges; j ++){
            if (edgesEqual(node1->edges[j], node2->edges[i]) == 1){
                edgeFound = 1;
                break;
            }

            // if the edge wasn't found, retur0
            if (edgeFound == 0) return 0;
        }
    }

    return 1;
}

/**
 * Recursively show this node and all of its descendants
 * 
 * @param mdn MassDawgNode *    the node to show and its descendants
 * @param spaces int            the number of spaces to print before
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
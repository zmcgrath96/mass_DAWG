typedef struct massDawgNode {
    struct Edge ** edges;
    int numEdges;
    char ** kmers;
    int final;
} MassDawgNode;

struct Edge {
    MassDawgNode * child;
    double mass;
};
/**
 * Allocate memory for a new MassDawgNode and init final to 0
 * 
 * @return MassDawgNode* to the new mass node
 * 
*/
MassDawgNode * initMassDawgNode() {
    MassDawgNode * mdn = (MassDawgNode *)malloc(sizeof(MassDawgNode));
    mdn->final = 0;
    mdn->numEdges = 0;
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
        int i;
        
        for (i = 0; i < mdn->numEdges; i++){
            deleteMassDawgNode(mdn->edges[i]->child);
            free(mdn->edges[i]);
        }
    }

    free(mdn);
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
    char ** newKmers = (char **)realloc(mdn->kmers, sizeof(mdn->kmers) + sizeof(char*));
    
    // add the new kmers to the end of the kmers
    int numKmers = sizeof(newKmers)/sizeof(char *);
    newKmers[numKmers-1] = kmer;

    // add it to the node
    mdn->kmers = newKmers;
}

/**
 * Add a new edge with a kmer and mass to the node
 * 
 * @param mdn MassDawgNode *         the node to add an edge to
 * @param kmer char *               the new kmer associated with the new node
 * @param mass double                the new mass to associate with the new node
 * 
 * @return None
*/
void addNewNode(MassDawgNode * mdn, char * kmer, double mass){

    // create the new mass node
    MassDawgNode * newNode = initMassDawgNode();
    newNode->kmers = malloc(sizeof(char*));
    newNode->kmers[0] = kmer;

    // make a new edge
    struct Edge * newEdge = (struct Edge *)malloc(sizeof(struct Edge *));
    newEdge->mass = mass;
    newEdge->child = newNode;

    // add the edge to the current node
    mdn->edges = (struct Edge *)realloc(mdn->edges, sizeof(struct Edge *) * (mdn->numEdges + 1));
    mdn->edges[mdn->numEdges] = newEdge;
    mdn->numEdges ++;
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
        if (mdn->edges[i]->mass == mass){
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
        if (mdn->edges[i]->mass == mass){
            return mdn->edges[i]->child;
        }
    }

    return initMassDawgNode();
}


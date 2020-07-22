#include <stdio.h>
#include <string.h>

#include <MassDawgNode.c>
#include <utils.c>

struct PreviousSequence {
    double * sequence;
    MassDawgNode ** nodes;
    int length; 
};

typedef struct massDawg {
    MassDawgNode * root;
    struct PreviousSequence * previousSequence;
    MassDawgNode ** uncheckedNodes;
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
    struct PreviousSequence * ps = (struct PreviousSequence *)malloc(sizeof(struct PreviousSequence *));
    double sequence[length];
    MassDawgNode * nodes[length];
    ps->sequence = sequence;
    ps->nodes = nodes;
    ps->length = length;
    
    return ps;
}

/**
 * Create a new MassDawg
 * 
 * @returns MassDawg *
*/
MassDawg * initMassDawg(){
    MassDawg * md = (MassDawg *)malloc(sizeof(MassDawg *));
    md->root = initMassDawgNode();
    md->previousSequence = initPreviousSequence(0);

    return md;
}

/**
 * Combine all possible nodes that can be combined from here down
 * 
 * @param md MassDawg * the graph we are reducing
 * @param downTo int the level of the graph down to which we should look to combine nodes
*/
void minimize(MassDawg * md, int downTo){
    return;
}

/**
 * Insert a sequence of masses into the Dawg. The input 
 * to the insertion (for many insertions) must be sorted.
 * 
 * @param md MassDawg * the MassDawg struct to insert into
 * @param seqeunce double [] the list of doubles that make up 
 * @param kmer char ** the sequence of amino acids associated with these masses
*/
void insert(MassDawg * md, double * sequence, char ** kmer){
    // if the previous sequence is not none (from init) and the 
    // new sequence is > the old sequence, then we need to error
    if (md->previousSequence->length > 0 && compareDoubleArrays(sequence, md->previousSequence->sequence)){
        printf("ERROR: Input sequence must be less than the last sequence.");
        return;
    }

    // Continue with adding the sequnce since the sequence is greater than the last

    // keep track of how long the common prefix is between the old sequence and
    // the new sequence
    int commonPrefixLength = 0;

    // minimum iterative length between the two sequences
    int iterLength = MIN(DOUBLE_ARR_LEN(sequence), DOUBLE_ARR_LEN(md->previousSequence->sequence));

    // go through each of the values in the two sequences and compare values until
    // we either run out of values or the two differ
    for (int i = 0; i < iterLength; i++) {
        
        // check to see if values at the two sequences at position i are the same.
        // if not, break 
        if (sequence[i] != md->previousSequence->sequence[i]) break;

        // for the last seuqnce, add the kmer up to this point to that node
        char newKmer[i+1];
        strncpy(newKmer, kmer, i);
        newKmer[i] = '\0'; 

        // check to see if any of the strings in this node (from previous sequence)
        // are this kmer. If no, add it
        MassDawgNode * currentNodeInPrevious = md->previousSequence->nodes[i];
        int numKmersInNode = STR_ARR_LEN(currentNodeInPrevious->kmers);

        // keep track if the kmer is found
        int kmerFound = 0;

        for (int kmerIdx = 0; kmerIdx < numKmersInNode; i++){
            // check to see if the two strings are equal
            if (strcmp(newKmer, currentNodeInPrevious->kmers[i]) != 0) {
                kmerFound = 1;
                break;
            }
        }

        // if the kmer wasn't found, add it to the node
        if (kmerFound == 1) addKmer(newKmer, currentNodeInPrevious);

        // we've made it this far so we increment
        commonPrefixLength ++;
    }

    // combine all nodes we can combine
    minimize(md, commonPrefixLength);

    // the current node we need to start at so that we can add the suffix without repeats
    MassDawgNode * currentNode;

    //add the suffix, starting from the correct node mid-way through the graph
    if (md->numUncheckedNodes == 0) currentNode = md->root;
    // get the next node (should be one) of the last unchecked node
    else currentNode = md->uncheckedNodes[md->numUncheckedNodes - 1]->edges[0]->child;

}

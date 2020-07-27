  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MassDawg.c"

int main(){
    // make a mass dawg and add a mass to it
    MassDawg * md = initMassDawg();

    double singly[4] = {100.1, 200.2, 300.3, 400.4};
    double doubly[4] = {200.2, 400.4, 600.6, 800.8};

    // printf("Making a new node and adding a kmer to it...\n");
    // MassDawgNode * newNode = initMassDawgNode();
    // addKmer("ABC", newNode);

    // printf("Adding a child to the node\n");
    // addNewNode(newNode, "ABCD", 100.1, 200.2);

    // printf("Deleting the nodes...\n");
    // deleteMassDawgNode(newNode);

    printf("Inserting a singly and doubly sequence with kmer ABC into DAWG\n");
    insert(md, singly, doubly, "ABCD", 4);

    // go through the graph at this point and print things out to see if it inserted
    printf("");
    showDawg(md);

    // clean up and remove the dawg
    deleteMassDawgNode(md->root);
    free(md);
}
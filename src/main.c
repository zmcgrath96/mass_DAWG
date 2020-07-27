  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MassDawg.c"

int main(){
    // make a mass dawg and add a mass to it
    MassDawg * md = initMassDawg();

    double singly[3] = {100.1, 200.2, 300.3};
    double doubly[3] = {200.2, 400.4, 600.6};

    printf("Inserting a singly and doubly sequence with kmer ABC into DAWG\n");
    insert(md, singly, doubly, "ABC", 3);

    // go through the graph at this point and print things out to see if it inserted
    printf("");
    showDawg(md);

    double singly2[4] = {100.1, 200.2, 300.3, 400.4};
    double doubly2[4] = {200.2, 400.4, 600.6, 800.8};

    printf("\nInserting a singly and doubly sequence with kmer ABCD into DAWG\n");
    insert(md, singly2, doubly2, "ABCD", 4);

    printf("");
    showDawg(md);

    // clean up and remove the dawg
    deleteMassDawgNode(md->root);
    free(md);
}
  
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

    double singly2[4] = {100.1, 200.2, 300.3, 400.4};
    double doubly2[4] = {200.2, 400.4, 600.6, 800.8};

    printf("\nInserting a singly and doubly sequence with kmer ABCD into DAWG\n");
    insert(md, singly2, doubly2, "ABCD", 4);

    printf("");
    showDawg(md);

    double singly3[4] = {100.1, 201.21, 302.32, 400.4};
    double doubly3[4] = {200.2, 402.42, 604.64, 800.8};
    printf("\nInserting a singly and doubly sequence with kmer AXYD into DAWG");
    insert(md, singly3, doubly3, "AXYD", 4);

    // finish the dawg
    printf("\nFinishing the dawg\n");
    finish(md);

    printf("");
    showDawg(md);

    // clean up and remove the dawg
    deleteMassDawg(md);
}
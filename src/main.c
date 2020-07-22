  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <MassDawgNode.c>
#include <utils.c>

int main(){
    MassDawgNode * md = (MassDawgNode *)malloc(sizeof(MassDawgNode *));
    printf("Hello\n");

    addNewNode(md, "ABC", 100.2);
    printf("%f    %s\n", md->edges[0]->mass, md->edges[0]->child->kmers[0]);

    int hasChild = hasChildWithValue(md, 100.2);
    printf("has child: %i\n", hasChild);

    double arr1[3] = {10.2, 11.6, 12.0};
    double arr2[4] = {3.12, 100.4, 1233.90, 9000.0};

    int lessThan = compareDoubleArrays(arr2, arr1);
    printf("%i result of comparison\n", lessThan);
}
#include <string.h>

#define MAX(a, b)           ((a < b) ?  (b) : (a))
#define MIN(a, b)           ((a > b) ? (b) : (a))

// end of string array
#ifndef EOSA
#define EOSA "end"
#endif

/**
 * Returns 1 if the first array of doubles is less than the second
 * otherwise return 0
 * 
 * @param array1 double * the first array to test
 * @param array2 double * the second array to test
 * 
 * @returns int 1 if array1 < array2 else 0
*/
int compareDoubleArrays(double * array1, double * array2, int array1Length, int array2Length) {

    // element wise compare
    int iterLength = MIN(array1Length, array2Length);

    for (int i = 0; i < iterLength; i ++){
        if (array1[i] > array2[i]) return 0;
    }
    return array1Length > array2Length ? 0 : 1;
}

/**
 * Determine if a value is the EOSA character
 * 
 * @param str1  char *  
 * 
 * @return 1 if is end 0 otherwise
*/
int isEOSA(char * str){
    if (str == NULL) return 1;
    return strcmp(str, EOSA) == 0 ? 1 : 0;
}

/**
 * Calculate the tolerance allowed in Daltons from a ppm number and a mass
 * 
 * @param mass      double  the mass we need to calculate the tolerance for
 * @param ppmTol    int     the tolerance in ppm to calculate
 * 
 * @return double           the tolerance in Daltons calculated
*/
double ppmToDal(double mass, int ppmTol){
    return ((double)ppmTol / 1000000.0) * mass;
}

/**
 * Combine two arrays of strings and return the concatenation of the
 * second onto the first
 * 
 * @param strArray1 char **     the first string array
 * @param strArray2 char **     the second string array
 * 
 * @return char **              the concatenated string arrays
*/
char ** concatStrArrays(char ** strArray1, char ** strArray2){
    // check to see if array 1 is null, if so return str array 2
    if (strArray1 == NULL) return strArray2;
    if (strArray2 == NULL) return strArray1;

    int strArray1Len = 0;

    // allocate the new concatenated array
    char ** concat = NULL;

    // iterate through each and allocate, copy, and free values
    while (1){
        // if we have hit the EOSA token, break
        if (isEOSA(strArray1[strArray1Len])) break;

        // get the length of the str at the index
        int sizeOfStr = strlen(strArray1[strArray1Len]);

        // allocate space for a string of the size of the one we are copying
        char * copiedStr = malloc(sizeof(char) * (sizeOfStr + 1));
        strcpy(copiedStr, strArray1[strArray1Len]);

        // set the end to NULL
        copiedStr[strArray1Len] = '\0';
        free(strArray1[strArray1Len]);
        strArray1[strArray1Len] = NULL;

        // make room in the concatenated array
        if (concat == NULL) concat = malloc(sizeof(char *));
        else concat = realloc(concat, sizeof(char *) * (strArray1Len + 1));

        concat[strArray1Len] = copiedStr;
        strArray1Len ++;
    }

    int i = 0;
    // copy things over from the last second array
    while (1){
        // if we hit the EOSA token, break
        if (isEOSA(strArray2[i])) break;

        // get the length of the string
        int sizeOfStr = strlen(strArray2[i]);
        
        // allocate space for a string the size of the one we are copying
        char * copiedStr = malloc(sizeof(char) * (sizeOfStr + 1));
        strcpy(copiedStr, strArray2[i]);

        // end the str with Null
        copiedStr[sizeOfStr] = '\0';
        free(strArray2[i]);
        strArray2[i] = NULL;

        // make room in the concatenated array
        // make room in the concatenated array
        if (concat == NULL) concat = malloc(sizeof(char *));
        else concat = realloc(concat, sizeof(char *) * (strArray1Len + i + 1));

        concat[strArray1Len + i] = copiedStr;
        i ++;
    }

    free(strArray1);
    free(strArray2);
    strArray1 = NULL;
    strArray2 = NULL;

    // make room for the EOSA token and add it to the end of the sequence
    concat = realloc(concat, sizeof(char *) * (strArray1Len + i + 2));
    concat[strArray1Len + i + 1] = EOSA;
    return concat;
}

/**
 * Deep copy the value of the two lists
 * 
 * @param array     char ** the list of strings to copy
 * 
 * @return char **          the copy of the list
*/
char ** deepCopyStrArray(char ** array){
    int arrayLen = (int)sizeof(array)/sizeof(char *);
    char ** copy = malloc(arrayLen * sizeof(char*));

    for (int i = 0; i < arrayLen; i ++){
        int lenOfStr = (int)sizeof(array[i])/sizeof(char);
        char * copiedStr = malloc(lenOfStr * sizeof(char));
        strcpy(copiedStr, array[i]);
        copiedStr[lenOfStr] = '\0';

        copy[i] = copiedStr;
    }
    return copy;
}
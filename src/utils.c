#define MAX(a, b)           ((a < b) ?  (b) : (a))
#define MIN(a, b)           ((a > b) ? (b) : (a))
#define STR_LEN(str)        ((int)sizeof(str)/sizeof(char))
#define DOUBLE_ARR_LEN(arr) ((int)sizeof(arr)/sizeof(double *))
#define STR_ARR_LEN(arr)    ((int)sizeof(arr)/sizeof(char *))

/**
 * Returns 1 if the first array of doubles is less than the second
 * otherwise return 0
 * 
 * @param array1 double * the first array to test
 * @param array2 double * the second array to test
 * 
 * @returns int 1 if array1 < array2 else 0
*/
int compareDoubleArrays(double * array1, double * array2) {
    // get the count of the two arrays
    int array1Length = DOUBLE_ARR_LEN(array1);
    int array2Length = DOUBLE_ARR_LEN(array2);

    // element wise compare
    int iterLength = MIN(array1Length, array2Length);

    for (int i = 0; i < iterLength; i ++){
        if (array1[i] > array2[i]) return 0;
    }
    return array1Length > array2Length ? 0 : 1;
}
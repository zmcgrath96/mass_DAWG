/**
 * Used for qsort on doubles. If d1 < d2, a number < 0 returned
 * if d1 > d2, a number > 0 returned
 * 
*/
int dblCmp(const void * d1, const void * d2){
    return *(double *)d1 - *(double *)d2;
}
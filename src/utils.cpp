/**
 * Used for qsort on doubles. If d1 < d2, a number < 0 returned
 * if d1 > d2, a number > 0 returned
 * 
*/
int dblCmp(const void * d1, const void * d2){
    return *(double *)d1 - *(double *)d2;
}

/**
 * Convert the ppm value of a given mass to a value in daltons
 * 
 * @param mass      double  the mass to calculate the dalton tolerance for
 * @param ppmTol    int     the tolerance in ppm
 * 
 * @return double   the tolerance in daltons 
*/
double ppmToDa(double mass, int ppmTol){
    return ((double)ppmTol / 1000000.0) * mass;
}
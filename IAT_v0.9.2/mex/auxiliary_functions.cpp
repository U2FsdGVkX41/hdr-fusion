
/*
 * This is a set of simple auxilliary functions towards basic operations
 * with matrices
 *
 * Copyright (C), 2012-2014 Georgios Evangelidis
 *
 * This file is part of the IAT library and is made available under
 * the terms of the GNU license (see the COPYING file).
 */

#include <matrix.h>
#include <mex.h>
#include <math.h>

int maxIndexOfArr(const mxArray* src){
    /*
     * It returns the largest-value's index of SRC array
     */
    
    double *a, maxi;
    int maxInd;
    a = mxGetPr(src);
    
    mwSize H = mxGetM(src);
    mwSize W = mxGetN(src);
    
    if ((W!=1) && (H!=1))
        mexErrMsgTxt("maxIndexOfArr: Input must be 1d array.");
    
    int i;
    
    maxInd=0;
    maxi=a[0];
    
    
    for (i=1; i<W*H; i++){
        if (a[i]>maxi){
            maxi = a[i];
            maxInd = i;
        }
    }
    
    return maxInd;
    
}

/////////////////////////////////////////////////////////////////////////

int maxIndexOfArr2D(const mxArray* src){
    /*
     * It returns the largest-value's index of SRC array
     */
    
    double *a, maxi;
    int maxInd;
    a = mxGetPr(src);
    
    mwSize H = mxGetM(src);
    mwSize W = mxGetN(src);
    
    int i;
    
    maxInd=0;
    maxi=a[0];
    
    if (mxIsNaN(maxi)){
        // check for nan
        for (i=1; i<W*H; i++){
            maxInd = i;
            maxi = a[i];
            if (!mxIsNaN(maxi)){
                break;
            }
        }
    }
    
    
    if (mxIsNaN(maxi)){
        mexErrMsgTxt("function maxIndexOfArr2D: All-Nan Array");
    }
    else {
        for (i=1; i<W*H; i++){
            if ((a[i]>maxi)){
                maxi = a[i];
                maxInd = i;
            }
        }
    }
    
    return maxInd;
    
}

/////////////////////////////////////////////////////////////////////////


int minIndexOfArr(const mxArray* src){
    /*
     * It returns the minimum-value's index of SRC array
     */
    
    double *a, maxi;
    int maxInd;
    a = mxGetPr(src);
    
    mwSize H = mxGetM(src);
    mwSize W = mxGetN(src);
    
    if ((W!=1) && (H!=1))
        mexErrMsgTxt("minIndexOfArr: Input must be 1d array.");
    
    int i;
    
    maxInd=0;
    maxi=a[0];
    
    
    for (i=1; i<W*H; i++){
        if (a[i]<maxi){
            maxi = a[i];
            maxInd = i;
        }
    }
    
    return maxInd;
    
}

/////////////////////////////////////////////////////////////////////////

double meanOfArr(const mxArray* src){
    /*
     * It returns the mean value of SRC array
     */
    
    double *a, average;
    a = mxGetPr(src);
    
    mwSize H = mxGetM(src);
    mwSize W = mxGetN(src);
    
    
    int i, j;
    
    average=0;
    
    for (i=0; i<W; i++){
        for (j=0; j<H; j++){
            average += a[i*H+j];
        }
    }
    average = average/(W*H);
    
    return average;
    
}

/////////////////////////////////////////////////////////////////////////


double normOfArr(const mxArray* src){
    /*
     * It returns the Euclidean norm of SRC array
     */
    
    double *a, norma;
    a = mxGetPr(src);
    
    mwSize H1 = mxGetM(src);
    mwSize W1 = mxGetN(src);
    
    
    int i, j;
    
    norma=0;
    
    for (i=0; i<W1; i++){
        for (j=0; j<H1; j++){
            norma +=a[i*H1+j]*a[i*H1+j];
        }
    }
    
    return sqrt(norma);
    
}


//////////////////////////////////////////////////////////////////////


void subArr(const mxArray* src1, const mxArray* src2, mxArray* dst){
    /*
     * It returns in DST array the difference SRC1-SRC2,
     * i.e. DST=SRC1-SRC2;
     */
    
    double *a, *b, *c;
    a = mxGetPr(src1);
    b = mxGetPr(src2);
    c = mxGetPr(dst);
    mwSize H = mxGetM(src1);
    mwSize W = mxGetN(src1);
    mwSize H2 = mxGetM(src2);
    mwSize W2 = mxGetN(src2);
    
    if ((W != W2) || (H != H2)) {
        mexErrMsgTxt("subArr: Input arrays must be of the same size.");
    }
    
    int i, j;
    
    for (i=0; i<W; i++){
        for (j=0; j<H; j++){
            c[i*H+j] = a[i*H+j]-b[i*H+j];
        }
    }
}

////////////////////////////////////////////////////////////


void addArr(const mxArray* src1, const mxArray* src2, mxArray* dst){
    /*
     * It returns in DST array the sum SRC1+SRC2
     */
    
    double *a, *b, *c;
    a = mxGetPr(src1);
    b = mxGetPr(src2);
    c = mxGetPr(dst);
    mwSize H = mxGetM(src1);
    mwSize W = mxGetN(src1);
    
    mwSize H2 = mxGetM(src2);
    mwSize W2 = mxGetN(src2);
    
    if ((W != W2) || (H != H2)) {
        mexErrMsgTxt("addArr: Input arrays must be of the same size.");
    }
    
    
    int i, j;
    
    for (i=0; i<W; i++){
        for (j=0; j<H; j++){
            c[i*H+j] = a[i*H+j]+b[i*H+j];
        }
    }
}




////////////////////////////////////////////////////////////


void addArrS(mxArray* src, const double offset){
    /*
     * It adds in SRC array the scalar OFFSET
     */
    
    double *a;
    a = mxGetPr(src);
    mwSize H = mxGetM(src);
    mwSize W = mxGetN(src);
    
    int i, j;
    
    for (i=0; i<W; i++){
        for (j=0; j<H; j++){
            a[i*H+j] += offset;
        }
    }
}


////////////////////////////////////////////////////////////


void mulArrS(mxArray* src, const double scale){
    /*
     * It mutliplies SRC array with scalar SCALE
     */
    
    double *a;
    a = mxGetPr(src);
    mwSize H = mxGetM(src);
    mwSize W = mxGetN(src);
    
    int i, j;
    
    for (i=0; i<W; i++){
        for (j=0; j<H; j++){
            a[i*H+j] *= scale;
        }
    }
}

////////////////////////////////////////////////////////////


void subArrS(mxArray* src, const double offset){
    /*
     * It subtracts the scalar OFFSET from SRC array
     */
    
    double *a;
    a = mxGetPr(src);
    mwSize H = mxGetM(src);
    mwSize W = mxGetN(src);
    
    int i, j;
    
    for (i=0; i<W; i++){
        for (j=0; j<H; j++){
            a[i*H+j] -= offset;
        }
    }
}


////////////////////////////////////////////////////////////

void normalizeMatrix(const mxArray* src, mxArray* dst){
    /*
     * It returns in DST array the SRC array normalized by its Euclidean
     * norm.
     */
    
    double *a, *b, norma;
    a = mxGetPr(src);
    b = mxGetPr(dst);
    
    mwSize H = mxGetM(src);
    mwSize W = mxGetN(src);
    int i, j;
    
    norma = 0.0;
    for (i=0; i<W; i++){
        for (j=0; j<H; j++){
            norma+= a[i*H+j]*a[i*H+j];
        }
    }
    
    norma = sqrt(norma);
    
    for (i=0; i<W; i++){
        for (j=0; j<H; j++){
            b[i*H+j]=a[i*H+j]/norma;
        }
    }
}

////////////////////////////////////////////////////////////

double dot2d(const mxArray* src1, const mxArray* src2){
    /*
     * It returns the dot product of vectors SRC1 and
     * SRC2 (when they are matrices the dot product is the one  of their
     * vectorized forms)
     */
    
    double *a, *b, out;
    a = mxGetPr(src1);
    b = mxGetPr(src2);
    
    mwSize H1 = mxGetM(src1);
    mwSize W1 = mxGetN(src1);
    
    mwSize H2 = mxGetM(src2);
    mwSize W2 = mxGetN(src2);
    
    if ((W1 != W2) || (H1 != H2)) {
        mexErrMsgTxt("dot2d:Inputs must be of the same size.");
    }
    
    int i, j;
    
    out=0;
    
    for (i=0; i<W1; i++){
        for (j=0; j<H1; j++){
            out +=a[i*H1+j]*b[i*H1+j];
        }
    }
    return out;
    
}

////////////////////////////////////////////////////////////

double sad(const mxArray* src1, const mxArray* src2){
    /*
     * It returns the sum-of-absolute differences of vectors SRC1 and
     * SRC2 (when they are matrices the dot product is that of their
     * vectorized forms)
     */
    
    double *a, *b, out;
    a = mxGetPr(src1);
    b = mxGetPr(src2);
    
    mwSize H1 = mxGetM(src1);
    mwSize W1 = mxGetN(src1);
    
    mwSize H2 = mxGetM(src2);
    mwSize W2 = mxGetN(src2);
    
    if ((W1 != W2) || (H1 != H2)) {
        mexErrMsgTxt("sad: Inputs must be of the same size.");
    }
    
    int i, j;
    
    out=0;
    
    for (i=0; i<W1; i++){
        for (j=0; j<H1; j++){
            out +=fabs(a[i*H1+j]-b[i*H1+j]);
        }
    }
    return out;
    
}

////////////////////////////////////////////////////////////

double ssd(const mxArray* src1, const mxArray* src2){
    /*
     * It returns the sum-of-squared-differences of vectors SRC1 and
     * SRC2 (when they are matrices the dot product is that of their
     * vectorized forms)
     */
    
    double *a, *b, out;
    a = mxGetPr(src1);
    b = mxGetPr(src2);
    
    mwSize H1 = mxGetM(src1);
    mwSize W1 = mxGetN(src1);
    
    mwSize H2 = mxGetM(src2);
    mwSize W2 = mxGetN(src2);
    
    if ((W1 != W2) || (H1 != H2)) {
        mexErrMsgTxt("ssd: Inputs must be of the same size.");
    }
    
    int i, j;
    
    out=0;
    double diff;
    
    for (i=0; i<W1; i++){
        for (j=0; j<H1; j++){
            diff = a[i*H1+j]-b[i*H1+j];
            out +=(diff*diff);
        }
    }
    return out;
    
}
////////////////////////////////////////////////////////////


void zeroPadding(const mxArray* src, mxArray* dst, int padW, int padH){
    /*
     * It returns in dst array the src array padded with zeros, such that
     * width(dst)=width(src)+2*padW and height(dst)=height(src)+2*padH
     */
    
    double *a, *b;
    a = mxGetPr(src);
    b = mxGetPr(dst);
    mwSize dimy = mxGetM(src);
    mwSize dimx = mxGetN(src);
    int i, j, fx, fy;
    
    fx = 2*padW+1;
    fy = 2*padH+1;
    
    
    for(i=0;i<dimx+fx-1;i++) {
        for(j=0;j<dimy+fy-1;j++) {
            if (i<padW || i>(dimx+padW-1) || j<padH || j>(dimy+padH-1)){
                b[i*(dimy+fy-1)+j] = 0;
            }
            else {
                b[i*(dimy+fy-1)+j] = a[(i-padW)*(dimy)+j-padH];
            }
        }
    }
    
    
}

////////////////////////////////////////////////////////////

void getPatch(const mxArray* src, mxArray* dst, int row, int col){
    /*
     * it returns in the DST array the patch of
     * the SRC array defined by the rectangular area
     * with origin (x,y)=(COL,ROW) and dimensions WIDTH and HEIGHT of DST
     */
    
    double *a, *b;
    a = mxGetPr(src);
    b = mxGetPr(dst);
    
    mwSize H1 = mxGetM(src);
    mwSize W1 = mxGetN(src);
    
    mwSize H2 = mxGetM(dst);
    mwSize W2 = mxGetN(dst);
    
    
    if ((row<0) || (col<0)) {
        mexErrMsgTxt("getPatch function: ROW and COL indices must be non-negative.");
    }
    
    if (((row+H2)>H1) || ((col+W2)>W1)) {
        mexErrMsgTxt("getPatch function: Patch extends outside image.");
    }
    
    int i, j;
    
    for(i=0;i<W2;i++) {
        for(j=0;j<H2;j++) {
            b[i*H2+j] = a[(col+i)*H1+(row)+j];
        }
    }
}

////////////////////////////////////////////////////////////

void getPatchZM(const mxArray* src, mxArray* dst, int row, int col){
    
    /* it returns in the DST array the zero-mean patch of
     * the SRC array (i.e. after subtracting its average)
     * defined by the rectangular area with origin
     * (x,y)=(COL,ROW) and dimensions WIDTH and HEIGHT of DST
     */
    
    double *a, *b;
    a = mxGetPr(src);
    b = mxGetPr(dst);
    
    mwSize H1 = mxGetM(src);
    mwSize W1 = mxGetN(src);
    
    mwSize H2 = mxGetM(dst);
    mwSize W2 = mxGetN(dst);
    
    
    if ((row<0) || (col<0)) {
        mexErrMsgTxt("getPatch function: ROW and COL indices must be non-negative.");
    }
    
    if (((row+H2)>H1) || ((col+W2)>W1)) {
        mexErrMsgTxt("getPatch function: Patch extends outside image.");
    }
    
    int i, j;
    double v = 0.0;
    
    for(i=0;i<W2;i++) {
        for(j=0;j<H2;j++) {
            
            b[i*H2+j] = a[(col+i)*H1+(row)+j];
            v += a[(col+i)*H1+(row)+j];
            
        }
    }
    
    v = v/(W2*H2);
    
    for(i=0;i<W2;i++) {
        for(j=0;j<H2;j++) {
            
            b[i*H2+j] = b[i*H2+j]-v;
            
        }
    }
    
    
}


///////////////////////////////////////////////////////////////////////////

void fillROI(mxArray* src, int row, int col, int width, int height, double entryValue) {
    /*
     * It fills a rectangular area of SRC array with value defined by ENTRYVALUE.
     * The rectangular area is defined by its top-left origin (ROW,COL) and
     * its WIDTH and HEIGHT.
     */
    
    double *a;
    a = mxGetPr(src);
    
    mwSize H = mxGetM(src);
    mwSize W = mxGetN(src);
    
    if ((H<(row+height)) || (W<(col+width)))
        mexErrMsgTxt("fillROI: the rectangular area extends out of array.");
    
    
    int i, j;
    
    for(i=0;i<width;i++) {
        for(j=0;j<height;j++) {
            a[(col+i)*H+(row)+j] = entryValue;
            //printf("%d\n",(col+i)*H+row+j);
        }
    }
}


///////////////////////////////////////////////////////////////////////////

void setZero(mxArray* src) {
    /*
     * It fills the SRC array with zeros.
     */
    
    double *a;
    a = mxGetPr(src);
    
    mwSize H = mxGetM(src);
    mwSize W = mxGetN(src);
    
    
    int i, j;
    
    for(i=0;i<W;i++) {
        for(j=0;j<H;j++) {
            a[i*H+j] = 0.0;
        }
    }
}

///////////////////////////////////////////////////////////////////////////

void setValue(mxArray* src, double entryValue) {
    /*
     * It fills the SRC with ENTRYVALUE.
     */
    
    double *a;
    a = mxGetPr(src);
    
    mwSize H = mxGetM(src);
    mwSize W = mxGetN(src);
    
    
    int i, j;
    
    for(i=0;i<W;i++) {
        for(j=0;j<H;j++) {
            a[i*H+j] = entryValue;
        }
    }
}
///////////////////////////////////////////////////////////////////////////

void setOne(mxArray* src) {
    /*
     * It fills the SRC array with ones.
     */
    double *a;
    a = mxGetPr(src);
    
    mwSize H = mxGetM(src);
    mwSize W = mxGetN(src);
    
    
    int i, j;
    
    for(i=0;i<W;i++) {
        for(j=0;j<H;j++) {
            a[i*H+j] = 1.0;
        }
    }
}

///////////////////////////////////////////////////////////////////////////

void setIdentity(mxArray* src) {
    /*
     * It tranlsates the square marix SRC into the identity matrix
     */
    
    double *a;
    a = mxGetPr(src);
    
    mwSize H = mxGetM(src);
    mwSize W = mxGetN(src);
    
    if (W!=H)
        mexErrMsgTxt("setIdentity: Array must be square.");
    
    int i, j;
    
    for(i=0;i<W;i++) {
        for(j=0;j<H;j++) {
            
            if (i==j){
                a[i*H+j] = 1.0;}
            else {
                a[i*H+j] = 0.0;}
        }
    }
}

///////////////////////////////////////////////////////////////////////////

void compArr(const mxArray* src1, const mxArray* src2, mxArray* dst) {
    /*
     * It compares SRC1 with SRC2 array and returns the binary array DST.
     * DST(i,j)=1 when SRC1(i,j)>SRC2(i,j), and 0 otherwise.
     */
    
    double *a, *b, *c;
    a = mxGetPr(src1);
    b = mxGetPr(src2);
    c = mxGetPr(dst);
    
    
    mwSize H1 = mxGetM(src1);
    mwSize W1 = mxGetN(src1);
    
    mwSize H2 = mxGetM(src2);
    mwSize W2 = mxGetN(src2);
    
    if ((W1 != W2) || (H1 != H2)) {
        mexErrMsgTxt("compArr: Compared array have different size.");
    }
    
    int i, j;
    
    for(i=0;i<W1;i++) {
        for(j=0;j<H1;j++) {
            
            if (a[i*H1+j]>b[i*H1+j]){
                c[i*H1+j] = 1.0;}
            else {
                c[i*H1+j] = 0.0;}
        }
    }
}

///////////////////////////////////////////////////////////////////////////

void compArrS(const mxArray* src1, double value, mxArray* dst) {
    /*
     * It compares elements of SRC1 with VALUE and returns the binary array DST.
     * DST(i,j)=1 when SRC1(i,j)>value, and 0 otherwise.
     */
    
    double *a, *c;
    a = mxGetPr(src1);
    c = mxGetPr(dst);
    
    
    mwSize H1 = mxGetM(src1);
    mwSize W1 = mxGetN(src1);
    
    mwSize H2 = mxGetM(dst);
    mwSize W2 = mxGetN(dst);
    
    if ((W1 != W2) || (H1 != H2)) {
        mexErrMsgTxt("compArrS: Source and target arrays have different size.");
    }
    
    int i, j;
    
    for(i=0;i<W1;i++) {
        for(j=0;j<H1;j++) {
            
            if (a[i*H1+j]>value){
                c[i*H1+j] = 1.0;}
            else {
                c[i*H1+j] = 0.0;}
        }
    }
}

///////////////////////////////////////////////////////////////////////////

void fillArr(mxArray* src, const mxArray* bin, double entryValue) {
    /*
     * It sets the values of SRC array indexed by the binary array BIN
     * equal to ENTRYVALUE
     */
    
    double *a, *b;
    a = mxGetPr(src);
    b = mxGetPr(bin);
    
    mwSize H1 = mxGetM(src);
    mwSize W1 = mxGetN(src);
    
    mwSize H2 = mxGetM(bin);
    mwSize W2 = mxGetN(bin);
    
    if ((W1 != W2) || (H1 != H2)) {
        mexErrMsgTxt("fillArr: Input and binary array have different size.");
    }
    
    int i, j;
    
    for(i=0;i<W1;i++) {
        for(j=0;j<H1;j++) {
            
            if (b[i*H1+j]==1.0)
                a[i*H1+j] = entryValue;
            
        }
    }
    
}


///////////////////////////////////////////////////////////////////////////

void equalArr(const mxArray* src, const mxArray* bin, mxArray* dst) {
    /*
     * It sets the values of DST array but indexed by the binary array BIN
     * equal to those of SRC array, i.e. DST(i,j)=SRC(i,j) if BIN(i,j)=1
     */
    
    double *a, *b, *c;
    a = mxGetPr(src);
    b = mxGetPr(bin);
    c = mxGetPr(dst);
    
    mwSize H1 = mxGetM(src);
    mwSize W1 = mxGetN(src);
    
    mwSize H2 = mxGetM(bin);
    mwSize W2 = mxGetN(bin);
    
    mwSize H3 = mxGetM(dst);
    mwSize W3 = mxGetN(dst);
    
    
    if ((W1 != W2) || (H1 != H2) || (H1!=H3) || (W1!=W3)) {
        mexErrMsgTxt("equalArr: Input arrays have different size.");
    }
    
    
    int i, j;
    
    for(i=0;i<W1;i++) {
        for(j=0;j<H1;j++) {
            
            if (b[i*H1+j]==1.0)
                c[i*H1+j] = a[i*H1+j];
            
        }
    }
    
}

////////////////////////////////////////////////////////////

void displayElements(const mxArray* src){
    /*
     * It displays the value of the SRC Array (first column first)
     */
    mexPrintf("1-based indexing\n");
    double *a;
    a = mxGetPr(src);
    
    mwSize H = mxGetM(src);
    mwSize W = mxGetN(src);
    int i, j;
    
    for (i=0; i<W; i++){
        for (j=0; j<H; j++){
            mexPrintf ("value(%d, %d): %f \n", j+1,i+1,a[i*H+j]);
        }
    }
}

////////////////////////////////////////////////////////////

void printMemoryArray(double* M, int rows, int cols) {
    /*
     * Prints ROWS*COLS elements starting from M's pointer as a Matlab
     * matrix of size ROWSxCOLS.
     */
    
    mxArray* matrix = mxCreateDoubleMatrix(rows, cols, mxREAL);
    double* matrixPtr = mxGetPr(matrix);
    int i;
    /* memcpy(matrixcells, M, rows*cols*sizeof(double)); */
    for (i = 0; i < rows*cols; i++) {
        *(matrixPtr++) = *(M++);
    }
    mexCallMATLAB(0, NULL, 1, &matrix, "disp");
}

/*
 * void printMexArray(const doublemxArray* M) {
 *
 * //  Prints the array M as in Matlab command line
 *
 * mexCallMATLAB(0, NULL, 1, &M, "disp");
 * }
 */



////////////////////////////////////////////////////////////

void displayMean(const mxArray* src){
    
    /*
     * It displays the average value of the SRC array
     */
    
    double *a;
    a = mxGetPr(src);
    
    mwSize H1 = mxGetM(src);
    mwSize W1 = mxGetN(src);
    
    int i, j;
    double v = 0.0;
    
    for(i=0;i<W1;i++) {
        for(j=0;j<H1;j++) {
            
            v += a[(i)*H1+j];
            
        }
    }
    
    v = v/(W1*H1);
    mexPrintf("mean value: %f\n", v);
    
}

////////////////////////////////////////////////////////////

int s2i(const int height, const int rowIndex, const int colIndex){
    /*
     * It converts zero-based 2D sub-indices into zero-based 1D column-wise index
     */
    
    return (colIndex*height+rowIndex);
}

////////////////////////////////////////////////////////////

#ifdef __cplusplus
void i2s(int height, int index, int &row, int &col) //C++
{
    // It converts zero-based 1D index into zero-based 2D sub-indices
    
    col = (int)(index)/height;
    
    row = index-(col*height);
    
}
#else
void i2s(int height, int index, int *row, int *col)//C
{
    // It converts zero-based 1D index into zero-based 2D sub-indices
    
    *col = (int)(index)/height;
    
    *row = index-((*col)*height);
    
}
#endif

///////////////////////////////////////////////////////////////////////////

void copyArr(const mxArray* src1, mxArray* src2) {
    /*
     * It really copies SRC1 to SRC2.
     */
    
    double *a, *b;
    a = mxGetPr(src1);
    b = mxGetPr(src2);
    
    mwSize H1 = mxGetM(src1);
    mwSize W1 = mxGetN(src1);
    
    mwSize H2 = mxGetM(src2);
    mwSize W2 = mxGetN(src2);
    
    if ((W1 != W2) || (H1 != H2)) {
        mexErrMsgTxt("Arrays have different size.");
    }
    
    int i;
    
    for(i=0;i<W1*H1;i++) {
        b[i] = a[i];
    }
    
}


///////////////////////////////////////////////////////////////////////////

int countValueInArray(const mxArray* src1, double value) {
    /*
     * It returns the number of appearances of VALUE in SRC1 array
     */
    
    double *a;
    a = mxGetPr(src1);
    
    
    mwSize H1 = mxGetM(src1);
    mwSize W1 = mxGetN(src1);
    
    
    int count;
    
    for(int i=0;i<W1*H1;i++) {
        
        if (a[i]==value){
            count++;
        }
    }
    
    return count;
}

///////////////////////////////////////////////////////////////////////////

void transposeArr(const mxArray* src, mxArray* dst){
    /*
     * It returns in DST array the SRC-transpose;
     */
    
    double *a, *c;
    a = mxGetPr(src);
    c = mxGetPr(dst);
    mwSize H = mxGetM(src);
    mwSize W = mxGetN(src);
    mwSize H2 = mxGetM(dst);
    mwSize W2 = mxGetN(dst);
    
    if ((W != H2) || (H != W2)) {
        mexErrMsgTxt("transposeArr: arrays with non compatible size for trasnposition.");
    }
    
    int i, j;
    
    for (i=0; i<W; i++){
        for (j=0; j<H; j++){
            c[j*W+i] = a[i*H+j];
        }
    }
}

////////////////////////////////////////////////////////////

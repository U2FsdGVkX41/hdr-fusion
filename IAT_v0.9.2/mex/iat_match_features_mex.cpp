
#include "auxiliary_functions.cpp"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    
    if ( (nrhs<2) || (nrhs>3) )
        mexErrMsgTxt("iat_match_descriptors: Use two or three input arguments.");
    
    if (nlhs!=4)
        mexErrMsgTxt("iat_match_descriptors: Use four output arguments.");
    
    if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]))
        mexErrMsgTxt("iat_match_descriptors: Use inputs with double type only.");
    
    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]))
        mexErrMsgTxt("iat_match_descriptors: Use inputs with double type only.");
    
    double ratio = 1.;
    
    if (nrhs ==3) {
        ratio = mxGetScalar(prhs[2]);
    }
    
    
    int m1 = mxGetM(prhs[0]);
    int n1 = mxGetN(prhs[0]);
    int m2 = mxGetM(prhs[1]);
    int n2 = mxGetN(prhs[1]);
    
    
    if (n1 != n2) {
        mexErrMsgTxt("iat_match_descriptors: descriptor matrices should have the same number of columns.");
    }
    
    if (ratio<0 || ratio>1) {
        mexErrMsgTxt("iat_match_descriptors: ratio should be in range (0,1].");
    }
    
    
    mxArray* res;
    double* pr;
    int i,j,k;
    
    res = mxCreateDoubleMatrix(m1,m2, mxREAL);
    pr =  mxGetPr(res);
    
    mxArray* vec1;
    vec1 = mxCreateDoubleMatrix(1,n1, mxREAL);
    mxArray* vec2;
    vec2 = mxCreateDoubleMatrix(1,n1, mxREAL);
    
    // compute pairwise angles from cosine
    for (i = 0; i < m1; i++) {
        getPatch(prhs[0], vec1, i, 0);
        double norm1 = normOfArr(vec1);
        for (j=0; j<m2; j++) {
            getPatch(prhs[1], vec2, j, 0);
            double norm2 = normOfArr(vec2);
            pr[i+j*m1]=acos(dot2d(vec1, vec2)/(norm1*norm2+0.0000001));
        }
    }
    
    
    mxArray* mapping = mxCreateDoubleMatrix(m1,1, mxREAL);
    double* mappingPr = mxGetPr(mapping);
    setZero(mapping);
    
    
    mxArray* temp = mxCreateDoubleMatrix(1,m2, mxREAL);
    double* tempPr =  mxGetPr(temp);
    
    mxArray* temp2 = mxCreateDoubleMatrix(m1,1, mxREAL);
    double* temp2Pr =  mxGetPr(temp2);
    
    mxArray* tempHist = mxCreateDoubleMatrix(1,m2, mxREAL);
    double* tempHistPr =  mxGetPr(tempHist);
    setZero(tempHist);
    
    
    //get potential matches
    int potmatch=0;
    for (i = 0; i < m1; i++) {
        getPatch(res, temp, i, 0);
        int index1 = minIndexOfArr(temp);
        double min1 = tempPr[index1];
        tempPr[index1] = 3.15; //pi+something
        int index2 = minIndexOfArr(temp);
        double min2 = tempPr[index2];
        
        if ( (min1<(ratio*min2)) && (pr[s2i(m1,i,index2)]>0)) {
            mappingPr[i] = index1; //zero-based indexing
            tempHistPr[index1] = tempHistPr[index1]+1;
        }
        
    }
    
    mxArray* pos1 = mxCreateDoubleMatrix(m1,1, mxREAL);
    double* pos1Pr =  mxGetPr(pos1);
    mxArray* pos2 = mxCreateDoubleMatrix(m1,1, mxREAL);
    double* pos2Pr =  mxGetPr(pos2);
    
    
    double value;
    for (i=0; i<m1;i++){
        value = tempHistPr[(int) mappingPr[i]];
        if (value==1){
            *pos1Pr = (double) i;
            *pos2Pr = (double) mappingPr[i];
            pos1Pr++;
            pos2Pr++;
            potmatch++;
        }
    }
    
    for (i=0; i<m2; i++){
        value = tempHistPr[i];
        if (value>1){
            
            for (j=0;j<m1;j++){
                if (mappingPr[j]==i)
                {
                    mappingPr[j] = 0;}
            }
            getPatch(res,temp2,0,i);
            int index1 = minIndexOfArr(temp2);
            if (mappingPr[index1] == 0){
                *pos1Pr = (double)index1;
                *pos2Pr = (double)i;
                pos1Pr++;
                pos2Pr++;
                potmatch++;
            }
        }
    }
    
    
    mexPrintf("Feature matching: %d correspondences found.\n",potmatch);
    if (potmatch<10)
        mexPrintf("iat_match_descriptors_mex: Increase the ratio to get more correspondences\n");
    
    
    pos1Pr -= potmatch;
    pos2Pr -= potmatch;
    
    plhs[2] = mxCreateDoubleMatrix(potmatch,1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(potmatch,1, mxREAL);
    
    double* pr2 = mxGetPr(plhs[2]);
    double* pr3 = mxGetPr(plhs[3]);
    
    plhs[0] = mxCreateNumericMatrix(m1, 1, mxDOUBLE_CLASS, mxREAL);
    double* pr0 = mxGetPr(plhs[0]);
    for (i=0; i<potmatch; i++){
        pr2[i] = pos1Pr[i];
        pr3[i] = pos2Pr[i];
        pr0[(mwSize) pr2[i] ] = pr3[i]+1; //add 1 for 1-based indexing
    }
    
    
    plhs[1] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
    double* pr1 = mxGetPr(plhs[1]);
    pr1[0] = (double) potmatch;
    
    
//convert to 1-based indexing
    addArrS(plhs[2],1);
    addArrS(plhs[3],1);
    
    
//There is no need to destroy. Matlab Memory Manager takes care of it:)
    
}

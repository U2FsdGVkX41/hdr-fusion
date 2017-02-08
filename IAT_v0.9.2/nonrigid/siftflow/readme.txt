This is the software package of ECCV paper:

C. Liu, J. Yuen, A. Torralba, J. Sivic and W. T. Freeman. SIFT flow: dense correspondence across different scenes. ECCV 2008.

Please cite this paper if you use our code for your research paper.

There is a big change compared to the original paper. We have a coarse-to-fine implementation of SIFT flow matching which runs much faster than the original algorithm presented in the paper.

Please go to "cpp" subfolder and follow README.txt to compile the cpp files.


Ce Liu
celiu@mit.edu
CSAIL MIT
Jan 2009


SIFT flow code is kindly provided from Ce Liu to IAT.

The compilation of SIFTflow through IAT, i.e. by using the function 'iat_mex'
compiles the file SIFTflow.cpp and saves the MEX-file iat_SIFTflow_mex.mex* 
in <iatRoot>/mex folder.
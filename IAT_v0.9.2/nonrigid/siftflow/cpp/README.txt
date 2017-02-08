NOTE: This is an updated README file from Georgios Evangelidis
compared to README file of original SIFTflow code provided by Ce Liu.

The *.cpp and *.h files in this folder are used to mex a matlab file for SIFT flow matching. 
Use the following command to compile in Matlab (tested on version 7.6 or later). You should
setup your mex environment first by 'mex -setup'.

The compilation should be fine through 'iat_mex' function. If you want to compile the file outside IAT try

>> mex mexDiscreteFlow.cpp BPFlow.cpp Stochastic.cpp


It has been tested in Windows x64, Linux x64 and Mac OS 10.8. 
Precompiled version for all platforms are provided by IAT.

------------------------- Important -------------------------

The following should be automatically addressed, if you got SIFTflow from IAT. 
However, if you still have problems, please follow Ce's advice and manually (un)comment
the line 
#define _LINUX_MAC
on the file "project.h",depending on your platform (uncomment the line in Windows).

-------------------------------------------------------------

Please contact Ce Liu (celiu@mit.edu) should you encounter any problems in compiling or find any bugs.
You can also contact Georgios Evangelidis (georgios@iatool.net), if you have problems with SIFTflow
because of its inclusion to IAT.

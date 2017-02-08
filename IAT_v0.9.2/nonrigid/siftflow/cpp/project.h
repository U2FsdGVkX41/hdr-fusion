#pragma once
#include "stdio.h"


//if this causes any problem, just define _LINUX_MAC only when you work on these platforms, and comment it if using visual studion in windows
#ifdef __APPLE__
#define _LINUX_MAC
#endif

#ifdef __linux__
#define _LINUX_MAC
#endif


template <class T>
void _Release1DBuffer(T* pBuffer)
{
	if(pBuffer!=NULL)
		delete []pBuffer;
	pBuffer=NULL;
}

template <class T>
void _Rlease2DBuffer(T** pBuffer,size_t nElements)
{
	for(size_t i=0;i<nElements;i++)
		delete [](pBuffer[i]);
	delete []pBuffer;
	pBuffer=NULL;
}


#define _MATLAB

#ifdef _MATLAB
#include "mex.h"
#endif


#ifdef _LINUX_MAC

template <class T1,class T2>
T1 __min(T1 a, T2 b)
{
  return (a>b)?b:a;

  printf("linux detected\n");

}

template <class T1,class T2>
T1 __max(T1 a, T2 b)
{
  return (a<b)?b:a;
}

#endif

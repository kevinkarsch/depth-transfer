// mex function to integrate annotations with belief propagation
#include "mex.h"
#include "project.h"
#include "Image.h"
#include <iostream>

using namespace std;

double* outputMatrix2(mxArray*& plhs,int n1,int n2)
{
	int dims[2];
	dims[0]=n1;
	dims[1]=n2;
	plhs=mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
	return (double *)mxGetData(plhs);
}

double* outputMatrix3(mxArray*& plhs,int n1,int n2,int n3)
{
	int dims[3];
	dims[0]=n1;
	dims[1]=n2;
	dims[2]=n3;
	plhs=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
	return (double *)mxGetData(plhs);
}

double* outputMatrix4(mxArray*& plhs,int n1,int n2,int n3,int n4)
{
	int dims[4];
	dims[0]=n1;
	dims[1]=n2;
	dims[2]=n3;
	dims[3]=n4;
	plhs=mxCreateNumericArray(4,dims,mxDOUBLE_CLASS,mxREAL);
	return (double *)mxGetData(plhs);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	DImage Im,Im2;
	IntImage vx,vy;
	
    //mexErrMsgTxt("error!\n");
    
	if(nrhs!=4)
		mexErrMsgTxt("Four parameters are needed for warping!\n");
	
	
	Im.LoadMatlabImage(prhs[0]);
	Im2.LoadMatlabImage(prhs[1]);
	vx.LoadMatlabImage(prhs[2],false);
	vy.LoadMatlabImage(prhs[3],false);
	
	if(!vx.matchDimension(vy))
		mexErrMsgTxt("The dimensions of horizontal and vertical flow fields must be identical!\n");
	
	if(Im.width()!=vx.width() || Im.height()!=vx.height())
		mexErrMsgTxt("The image and the flow field must have identical width and height!\n");
		
	DImage warpIm(Im.width(),Im.height(),Im.nchannels());
	IntImage mask(Im.width(),Im.height());
	
	// the main function is here
	size_t width=Im.width(),height=Im.height(),nchannels=Im.nchannels();
	size_t width2 = Im2.width(), height2 = Im2.height();
	for(size_t i=0;i<height;i++)
		for(size_t j=0;j<width;j++)
		{
			size_t offset=(i*width+j)*nchannels;
			size_t x=j+vx.data()[i*width+j];
			size_t y=i+vy.data()[i*width+j];
			//cout<<"x: "<<x<<" y:"<<y<<endl;
			//mexPrintf("vx: %d  vy: %d\n",vx.data()[i*width+j],vy.data()
			//mexPrintf("x: %d  y: %d\n",x,y);
			if(x>=0 && x<=width2-1 && y>=0 && y<=height2-1)
			{
				size_t offset2=(y*width2+x)*nchannels;
				//for(size_t k=0;k<nchannels;k++)
				//	warpIm.data()[offset+k]=Im.data()[offset1+k];
				memcpy(warpIm.data()+offset,Im2.data()+offset2,sizeof(double)*nchannels);
				mask.data()[i*width+j]=1;
			}
		}
	
	// output to Matlab
	warpIm.OutputToMatlab(plhs[0]);
	mask.OutputToMatlab(plhs[1]);
}
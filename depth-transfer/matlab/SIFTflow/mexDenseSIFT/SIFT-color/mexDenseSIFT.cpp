#include "mex.h"
#include "Image.h"
#include "ImageFeature.h"
#include "Vector.h"
#include <vector>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	DImage im;
	im.LoadMatlabImage(prhs[0]);

	bool IsMultiScale = false;
	vector<int> cellSizeVect;
	int cellSize = 3;
	int stepSize = 1;
	double normalization = 0.1;
	bool IsBoundaryIncluded = true;
    double alpha = 1;
    double beta = 1;
    int winSize = 3;

	Vector para;
    // load alpha and beta
    if(nrhs>1)
    {
        para.readVector(prhs[1]);
        alpha = para[0];
        if(para.dim()>1)
            beta = para[1];
    }
    // load cell size for dense SIFT image
	if(nrhs>2) // if cell size is input
	{
		para.readVector(prhs[2]);
		if(para.dim()>1)
		{
			IsMultiScale = true;
			for(int i = 0;i<para.dim();i++)
				cellSizeVect.push_back(para[i]);
		}
		else
			cellSize = para[0];
	}
	if(nrhs>3)
	{
		para.readVector(prhs[3]);
		winSize = para[0];
        if(para.dim()>1)
            stepSize = para[1];
        if(para.dim()>2)
            normalization = para[2];
        if(para.dim()>3)
            IsBoundaryIncluded = para[3];
    }
    
    // compute the SIFT image
	UCImage imsift;
	if(IsMultiScale)
		ImageFeature::imSIFT(im,imsift,cellSizeVect,stepSize,normalization,IsBoundaryIncluded);
	else
		ImageFeature::imSIFT(im,imsift,cellSize,stepSize,normalization,IsBoundaryIncluded);
    if(beta == 0)
    {
        imsift.OutputToMatlab(plhs[0]);
        return;
    }
    // now compute the patch image
    im.Multiplywith(255);
    UCImage Im,impatch;
    Im.copy(im);
    ImageFeature::imPatch(Im,impatch,winSize);
    
    // merge these two images
    Im.Merge(imsift,impatch,alpha,beta);
    
	Im.OutputToMatlab(plhs[0]);
}
# Depth Transfer
This MATLAB software accompanies our TPAMI'14 and ECCV'12 papers:

K. Karsch, C. Liu and S.B. Kang. "DepthTransfer: Depth Extraction from Video Using Non-parametric Sampling." IEEE TPAMI, 2014.

K. Karsch, C. Liu and S.B. Kang. "Depth Extraction from Video Using 
Non-parametric Sampling." ECCV 2012.

Descriptions of parameters and implementation details can be found in the paper and supplementary material (available at [https://kevinkarsch.com](https://kevinkarsch.com)). Please cite the TPAMI paper 
if you use our code for your research.

## A note on the (defunct) MSR-V3D Dataset
Microsoft inadvertently removed several of their datasets and research pages during a recent server change, including MSR-V3D, and I have been told that it is unlikely ever be replaced. However, there are other datasets which are all-in-all better, serve similar purposes:

- [NYU-V2](http://cs.nyu.edu/~silberman/datasets/nyu_depth_v2.html) 
- [Make3D](http://make3d.cs.cornell.edu/data.html)
- [Sintel](http://sintel.is.tue.mpg.de/depth)
- [Kitti](http://www.cvlibs.net/datasets/kitti)


## Installation
Pre-mexed binaries for Windows (32b and 64b), OS X (64b), and Linux (32b 
and 64b) are packaged with the code. If you're using any of those OSes, 
you should not need to compile anything. These binaries and code have
been verified with MATLAB 7.14 (2012a). 

If you need to compile the code (either on a different OS or the binary 
isn't suitable for your version of MATLAB), all you need to do is run 
'make' from a MATLAB prompt. If you run into any errors, see make.m. 

IMPORTANT: to use optical flow (motion estimation + temporal consistency), 
you must also download [Ce Liu's publicly available optical flow code](http://people.csail.mit.edu/celiu/OpticalFlow/). Modify the path in opticalflow.m (`OPTICAL_FLOW_DIR`) to match the location 
that you unpacked the OpticalFlow directory (the default is to unpack to 
this directory). Note that if you are only running on single images or on 
data from our dataset (StereoRGBD Dataset), this is not necessary since the
flow has been precomputed.

## Getting started
To get started, take a look at demoDepthTransfer.m, or just run 
'demoDepthTransfer' from a MATLAB prompt. This runs depth transfer on a 
single image from the Make3D Range Image Dataset.

To evaluate our code on the Make3D dataset, first read the documentation 
in evaluateMake3dDataset.m and run 'evaluateMake3dDataset' from a MATLAB 
prompt.

For more advanced usage, there are several scripts in the examples 
directory (these require the optical flow module; see Installation above):

```
  example1.m - Demonstrates depth transfer on a static video sequence with 
               moving objects
  example2.m - Demonstrates depth transfer on a rotating video sequence
  example3.m - Demonstrates depth transfer on a sequence with varying  
               focal length and moving objects
  example4.m - If you have also downloaded our StereoRGBD dataset, this 
               example will run one of the results using our training/test
               split.
```

## Directory contents
```
DepthTransfer/            - Root directory
  README.txt              - this file
  make.m                  - Helper function to compile SIFTflow and other 
                            cpp-mex files
  createData.m            - Stores raw data in DepthTransfer format
  loadData.m              - Reads DepthTransfer data stored on disk
  initializeProject.m     - Sets up a DepthTransfer project (params, paths,
                            etc)
  depthTranser.m          - main DepthTransfer module (relies on code in 
                            matlab/)
  demoDepthTransfer.m     - Simple demonstration of the code on a single image
  evaluateMake3dDataset.m - reproduces Make3D results from our paper
  opticalflow.m           - A wrapper for Liu's optical flow (see above)
  fillDepthHoles.m        - Code to spatiotemporally interpolate holes in 
                            depth maps
  demo/                   - Scripts and data for running advanced 
                            DepthTransfer demos
  matlab/                 - Directory containing core DepthTransfer code
    Features/               - Code for computing image and flow descriptors
	Inference/              - Code for inferring depth from warped 
                              candidate depth images
	MotionEstimation/       - Several algorithms for estimating a binary 
                              motion segmentation in videos
	SIFTflow/               - A modified version of Liu et al.'s SIFTflow**
    Utilties/               - DepthTransfer helper functions
  examples/                 - Example scripts and data
  data/                     - Default data directory (empty initially)
  results/                  - Default results directory (empty initially)

**Original SIFTflow code and paper can be found here:
    http://people.csail.mit.edu/celiu/SIFTflow/
```
	
## Notes
This code is not optimized for efficiency in any sense, and can take both a 
  very long time to run and a great deal of memory. Take advantage of 
  MATLAB's parallel toolbox if available (ie. 'matlabpool'), and be careful
  about the size of the input sequence. Using a project size of height=300
  and width=200, a video of 60 frames can take several GB of RAM.
  

##  Bug reports
For bug reports, please send me an email at kevin@kevinkarsch.com

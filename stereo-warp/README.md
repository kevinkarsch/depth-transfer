This MATLAB software accompanies our TPAMI'14 and ECCV'12 papers:

K. Karsch, C. Liu and S.B. Kang. "DepthTransfer: Depth Extraction from 
Video Using Non-parametric Sampling." IEEE TPAMI, 2014.

K. Karsch, C. Liu and S.B. Kang. "Depth Extraction from Video Using 
Non-parametric Sampling." ECCV 2012.

Disclaimer: this code is not well documented and is intended for use 
within our DepthTransfer system to demonstrate stereo conversion!

Descriptions of parameters and implementation details can be found in the 
papers and supplementary material (available at 
http://kevinkarsch.com/depthtransfer). Please cite one of these papers if 
you use our code for your research.

-------

This code is for converting monocular RGBD images/frames into stereo RGB 
frames (e.g. for 3D/stereoscopic viewing).

Run demo_warp.m for a simple demo, and see stereoWarp.m for details. 
Running demo_warp will produce output in the results/ folder. Initially
the results folder contains the results the code should produce. The demo
uses estimated depth, so the results might not look perfect; it should 
look better with ground truth depth maps.

You can leave the parameter settings as they are, but if you'd like a 
smoother result, you can adjust the weights near the top of stereoWarp.m. 
To get more/less disparity in the stereo images, try changing the 
displevels parameter (higher => more pronounced 3D effect, but possibly 
more eye strain).


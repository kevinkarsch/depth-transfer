This directory contains three different motion segmentation implementations 
(intended for use with DepthTransfer code). The main files are:

segmentMotionSimple.m           - Estimates motion by thresholding optical 
                                  flow. Works best on short sequences (a 
                                  few frames) when the camera is static.
segmentMotionHeuristic.m        - Estimates motion using the heuristic 
                                  described in our ECCV'12 paper. Works
                                  best on longer videos shot with a static 
                                  camera.
segmentGlobalMotionSimple.m*    - Same as segmentMotionSimple, but first
                                  computes a background mosaic to detect
                                  only motion relative to the camera (i.e. 
                                  not changes in camera rotation or focal 
                                  length). Use this for sequences that have 
                                  a rotating of zooming camera.
segmentGlobalMotionHeuristic.m* - Same as segmentMotionHeuristic, but first
                                  computes a background mosaic to detect
                                  only motion relative to the camera (i.e. 
                                  not changes in camera rotation or focal 
                                  length). Works best for longer sequences 
                                  that have a rotating of zooming camera.


These are fairly basic implementations, and may break for complex 
sequences. Note that NONE of these functions will work well if the camera 
contains a non-trivial amount of translation. DepthTransfer code allows for 
other motion segmentation functions if you prefer a different algorithm 
(see depthTransfer.m). For implementation details, please refer to our 
ECCV'12 paper: 

K. Karsch, C. Liu and S.B. Kang. "Depth Extraction from Video Using 
Non-parametric Sampling." ECCV 2012.

*NOTE: segmentGlobalMotion{*} are heavily dependant on 3rd party code:
1) An optical flow module is required, and uses the opticalflow.m wrapper 
   (found in the root of the depthTransfer code). The wrapper provides
   details on downloading and installing Ce Liu's publicly available
   module, found here: http://people.csail.mit.edu/celiu/OpticalFlow/
2) A subset of Peter Kovesi's MATLAB functions for computer vision
   (http://www.csse.uwa.edu.au/~pk/research/matlabfns/) have been included
   in the directory 'kovesi_fns'. Please see Peter's website for the full
   and up to date version of his code.
3) A subset of Andrea Vedaldi's SIFT code and binaries for MATLAB 
   (http://www.vlfeat.org/~vedaldi/code/sift.html) has been included
   in the directory 'sift'. Please see Andrea's website for the full and up 
   to date version of his code.

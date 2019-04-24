This directory contains example scripts for running Depth Transfer code. 
The examples are as follows:

```
example1 - Depth Transfer for a static video containing moving objects
example2 - Depth Transfer for a video captured with a rotating camera 
example3 - Depth Transfer for a video with variable focal length (zoom)
           and moving objects
example4 - Example usage of StereoRGBD dataset (downloaded separately*)
           with Depth Transfer
           
```

The other directories contain sample training and testing data used in the
example scripts. NOTE: `sample_training_data/` and `demo_data/`contain several
files from the Make3D Laser+Image dataset. The full data can be downloaded
here: http://make3d.cs.cornell.edu/data.html.

Note that these examples use an extremely small training set (only 7 RGBD
pairs in `sample_training_data`). More training data created using our own
StereoRGBD dataset, or any existing RGBD dataset.

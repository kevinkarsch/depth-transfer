function flow = opticalflow(img1, img2, para)
%OPTICALFLOW Compues optical flow from img1 to img2. 
%
% Input:
%  img1 - [height x width x {1,3}] image
%  img2 - [height x width x {1,3}] image
%  para - Parameter vector (default defined below)
%
% Output:
%  flow - [height x width x 2] array of flow vectors. flow(:,:,1) encodes
%         motion in the horizontal direction, flow(:,:,2) encodes motion in
%         the vertical direction.
%
% NOTE:
% This code depends on publicly available optical flow code from Ce Liu, 
% found here: http://people.csail.mit.edu/celiu/OpticalFlow/.
%
% Download the code from the above link and set the below variable to the
% path that the OpticalFlow directory exists (default is within the same
% directory as opticaflow.m). Compile the code by following the readme 
% of Liu's optical flow package.
opticalflow_m_dir = fileparts(mfilename('fullpath')); %Absolute path of this function (opticalflow.m)
OPTICAL_FLOW_DIR = fullfile(opticalflow_m_dir,'OpticalFlow');
%
%%%%%%%%%%%   Begin opticalflow   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(nargin<3) %See Coarse2FineTwoFrames for parameter details
        alpha = 0.012; %Regularization
        ratio = 0.75; %Downsample ratio
        minWidth = 30; %Width of coarsest level
        nOuterFPIterations = 7;
        nInnerFPIterations = 1;
        nSORIterations = 30;
        para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];
    end
    addpath(fullfile(OPTICAL_FLOW_DIR, 'mex'));
    [vx,vy] = Coarse2FineTwoFrames(img1, img2, para);
    flow = cat(3,vx,vy);
end


function project = initializeProject(Cv, Cf, dims)
%INITIALIZEPROJECT Sets globals for a Depth Transfer project
%
% Input:
%  Cv(=7)            - Number of top candidate videos to use for inference
%                      during Depth Transfer inference
%  Cf(=1)            - Number of top candidate frames from each candidate 
%                      video to use for inference during Depth Transfer 
%                      inference. Thus, there are a total of Cv*Cf
%                      candidates used during inference (per frame of
%                      input)
%  dims(=[290,215])  - Spatial dimensions of all Depth Transfer results.
%                      Before candidate matching and inference, all 
%                      testing/training data will be scaled to match these
%                      dimensions
%
% Output:
%  project - Struct containing information about a Depth Transfer project
%
% NOTE: to change directory structure, modify the following
ROOT_DIR = fileparts(mfilename('fullpath')); %Set root to the same directory as initializeProject.m
DATA_DIR = fullfile(ROOT_DIR, 'data');
RESULT_DIR = fullfile(ROOT_DIR, 'results');
CODE_DIR = fullfile(ROOT_DIR, 'matlab'); %Make sure to move code if modified
%
%%%%%%%%%%%   Begin initializeProject   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Project paths
    project.path.data       = DATA_DIR; 
    project.path.results    = RESULT_DIR;
    if( ~exist(fullfile(project.path.results),'dir') )
        mkdir(fullfile(project.path.results));
    end
    addpath(genpath(CODE_DIR)); %Add code paths
    
    %Project size and kNN info
    try
        project.Cv = Cv;
        project.Cf = Cf;
        project.h = dims(1);
        project.w = dims(2);
    catch %#ok<CTCH>
        warning('Not enough parameters, using defaults.\n'); %#ok<WNTAG>
        project.Cv      = 7; %Number of candidate videos
        project.Cf      = 1; %Number of candidate frames (per candidate video)
        %Default height and width is half of the standard StereoRGBD image dimensions (579x430)
        project.h       = 290; %Height of all images (input and training data will be scaled)
        project.w       = 215; %Width of all images (input and training data will be scaled)
    end
    
    %Depth Transfer optimization parameters
    %Coefficient on transferred relative depth matching (as opposed to absolute depth matching)
    project.para.depthTransferParams.gradCoeff       = 10;
    %Spatial smoothness coefficent
    project.para.depthTransferParams.smoothCoeff_spa = 10;
    %Temporal consistency coefficent
    project.para.depthTransferParams.smoothCoeff_tem = 100;
    %Coeffiecent weights the influence that the prior has on inferred depth
    project.para.depthTransferParams.priorCoeff      = 0.5;
    %Coeffiecent weights the influence that the motion segmentation has on inferred depth
    project.para.depthTransferParams.motionCoeff     = 5;
	%Weights sift and flow matching error confidence weights (higher => more emphasis on flow error)
    project.para.depthTransferParams.alpha           = 1;
    %Mixing weight between flow and appearance for choosing candidate matching frames
    project.para.depthTransferParams.gf_mix          = 0.5; 
    
    %Dense SIFT parameters
    project.para.SIFT.cellsize       = 4;
    project.para.SIFT.gridspacing    = 1;
    
    %SIFT flow parameters
    project.para.SIFTflowpara.alpha          = 1*255;   
    project.para.SIFTflowpara.d              = project.para.SIFTflowpara.alpha*20*255;
    project.para.SIFTflowpara.gamma          = 0.001*255;
    project.para.SIFTflowpara.nlevels        = 4;
    project.para.SIFTflowpara.topwsize.x     = 20;
    project.para.SIFTflowpara.topwsize.y     = 5;
    project.para.SIFTflowpara.wsize.x        = 4;
    project.para.SIFTflowpara.wsize.y        = 2;
    project.para.SIFTflowpara.nIterations    = 80;
    project.para.SIFTflowpara.nTopIterations = 100;
end

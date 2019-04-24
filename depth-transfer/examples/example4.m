%EXAMPLE4 This example shows how generate results on the Indoor portion of
%         our StereoRGBD dataset downloaded separately at:
%                http://kevinkarsch.com/depthtransfer
%         This script requires either the full or Indoor part of the
%         dataset (RAW or compressed version is fine).
% 
% Before running, make sure the StereoRGBD dataset (outdoor portion not 
% required) has been downloaded and unpacked to the data directory 
% (default is [DEPTH_TRANSFER_ROOT]/data). The directory should contain: 
%    [DEPTH_TRANSFER_ROOT]/data/Indoor1A-001/   (First indoor data)
%                       ...
%    [DEPTH_TRANSFER_ROOT]/data/Indoor4-009/    (Last indoor data)
% To change the name/location of the data directory, see variable DATA_DIR 
% in initializeProject.m.
%
% NOTE: Running this script can take up to 8GB of RAM! Computing the prior
% and features will also take a long time on the first run, but these are
% stored on disk and will not need to be computed on future runs.
%
%%%%%%%%%%%   Begin example4   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..'); %Add the root of depthTransfer to the path

%% Initialize a new project (sets up paths, parameters, etc)
h = 290; w = 215; %Inferred depth resolution (output)
Cv = 7; %Number of candidate videos to use for training
Cf = 1; %Number of candidate frames from each video
project = initializeProject(Cv, Cf, [h,w]);

%% Depth transfer
%Train with all data from Building 1
trainFiles = dir(fullfile(project.path.data, 'Indoor1*'));
%%% Examples of creating other trainFile structs:
% %Train with data from Building 1 and 2 (example)
% trainFiles = [dir(fullfile(project.path.data, 'Indoor1*')); ...
%               dir(fullfile(project.path.data, 'Indoor2*'))];
% %Train with data from Building 1A, 1C, and 3 (example)
% trainFiles = [dir(fullfile(project.path.data, 'Indoor1A*')); ...
%               dir(fullfile(project.path.data, 'Indoor1C*')); ...
%               dir(fullfile(project.path.data, 'Indoor3*'))];

%For this example, we infer depth for the last data in Building 4, clip 001
testFile = fullfile('Indoor4-009', '001');

%Compute prior (average training depth). To save time for future runs, here
% we precompute and store the prior (training data must stay constant).
if( exist('Indoor1_prior.mat', 'file') )
    load('Indoor1_prior.mat');
else
    fprintf('Computing depth prior...'); testTime = tic;
    depthPrior = computePrior(project, trainFiles);
    save('Indoor1_prior.mat', 'depthPrior');
    fprintf('done.   [%6.02fs]\n', toc(testTime));
end

%Set motion segmentation function
motionFunc = @segmentMotionHeuristic; %This motion estimator works best on 
                                      % longer sequences with a static camera

%Run depth transfer
depthEst = depthTransfer(project, testFile, trainFiles, depthPrior, motionFunc);

%% Display results
imgs = loadData(fullfile(project.path.data,'Indoor4-009'),'001', [project.h, project.w]);
NdepthEst = repmat(imnormalize(depthEst),[1,1,3,1]); %Normalize/add channels for visualization
figure('Name','Input images and estimated depth'),
imshow([imgs(:,:,:,10), imgs(:,:,:,30), imgs(:,:,:,50); ...
        NdepthEst(:,:,:,10), NdepthEst(:,:,:,30), NdepthEst(:,:,:,50)]);
implay([imgs, NdepthEst]); %Play as video

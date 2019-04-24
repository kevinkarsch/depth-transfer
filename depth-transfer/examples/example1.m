%EXAMPLE1 This example shows how to use depth transfer on an input sequence
%         with a static viewpoint and moving objects
%
%%%%%%%%%%%   Begin example1   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..'); %Add the root of depthTransfer to the path

%% Initialize a new project (sets up paths, parameters, etc)
h = 290; w = 215; %Inferred depth resolution (output)
Cv = 7; %Number of candidate videos to use for training
Cf = 1; %Number of candidate frames from each video
project = initializeProject(Cv, Cf, [h,w]);

%% Create training data
fprintf('Preparing training data...'); prepTrainTime = tic;
add_training_data(project); %Look here for details on creating traing data
fprintf('done. [%6.02fs]\n', toc(prepTrainTime));

%% Create test data for example1
example1Files = dir( fullfile('example1_data', '*.jpg') ); %Get listing of all images
imgExample1Cell = cell(numel(example1Files),1); %Initialize image cell
for i=1:numel(example1Files) %Read all images as doubles
    imgExample1Cell{i} = im2double(imread(fullfile('example1_data',example1Files(i).name)));
end
imgExample1 = cell2mat(reshape(imgExample1Cell,1,1,1,[])); %Convert cell -> 4D array
%Add data (unless it already exists from a previous run)
if( ~exist(fullfile(project.path.data,'example1','001'), 'dir') )
    createData(fullfile(project.path.data,'example1'), imgExample1); %Reformat data
end

%% Depth transfer
%Set which data to train with
trainFiles = dir(fullfile(project.path.data, 'Make3D-Train*'));
%Infer depth for created data (data/example1/001)
testFile = fullfile('example1', '001');
%Compute prior (average training depth). To save time for future runs, here
% we precompute and store the prior (training data must stay constant).
if( exist('sample_training_prior.mat', 'file') )
    load('sample_training_prior.mat');
else
    fprintf('Computing depth prior...'); testTime = tic;
    depthPrior = computePrior(project, trainFiles);
    save('sample_training_prior.mat', 'depthPrior');
    fprintf('done.   [%6.02fs]\n', toc(testTime));
end
%Set motion segmentation function
motionFunc = @segmentMotionSimple; %This one works best on short sequences 
                                   % with a static viewpoint (no rotation,
                                   % translation, zooming, etc)
%Run depth transfer
depthEst = depthTransfer(project, testFile, trainFiles, depthPrior, motionFunc);

%% Display results
imgs = loadData(fullfile(project.path.data,'example1'),'001', [project.h,project.w]);
NdepthEst = repmat(imnormalize(depthEst),[1,1,3,1]); %Normalize/add channels for visualization
figure('Name','Input images and estimated depth'),
imshow([imgs(:,:,:,1), imgs(:,:,:,2), imgs(:,:,:,3); ...
        NdepthEst(:,:,:,1), NdepthEst(:,:,:,2), NdepthEst(:,:,:,3)]);
implay([imgs, NdepthEst]); %Play as video

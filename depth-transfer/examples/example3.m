%EXAMPLE3 This example shows how to use depth transfer on an input sequence
%         with a zooming camera (variable focal length) and moving objects
%
% NOTE: Running this script can take up to 8GB of RAM!
%
%%%%%%%%%%%   Begin example3   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..'); %Add the root of depthTransfer to the path

%% Initialize a new project (sets up paths, parameters, etc)
h = 240; w = 320; %Inferred depth resolution (output)
Cv = 7; %Number of candidate videos to use for training
Cf = 1; %Number of candidate frames from each video
project = initializeProject(Cv, Cf, [h,w]);

%% Create training data
fprintf('Preparing training data...'); prepTrainTime = tic;
add_training_data(project); %Look here for details on creating traing data
fprintf('done. [%6.02fs]\n', toc(prepTrainTime));

%% Create test data for example3
example3Files = dir( fullfile( 'example3_data', '*.jpg') ); %Get listing of all images
imgExample3Cell = cell(numel(example3Files),1); %Initialize image cell
for i=1:numel(example3Files) %Read all images as doubles
    imgExample3Cell{i} = im2double(imread(fullfile('example3_data',example3Files(i).name)));
end
imgExample3 = cell2mat(reshape(imgExample3Cell,1,1,1,[])); %Convert cell -> 4D array
%Add data (unless it already exists from a previous run)
if( ~exist(fullfile(project.path.data,'example3','001'), 'dir') )
    createData(fullfile(project.path.data,'example3'), imgExample3); %Reformat data
end

%% Depth transfer
%Set which data to train with
trainFiles = dir(fullfile(project.path.data, 'Make3D-Train*'));
%Infer depth for created data (data/example3/001)
testFile = fullfile('example3', '001');
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
motionFunc = @segmentGlobalMotionHeuristic; %This one works best on longer sequences 
                                             % with a rotating/zooming camera
%Run depth transfer
depthEst = depthTransfer(project, testFile, trainFiles, depthPrior, motionFunc);

%% Display results
imgs = loadData(fullfile(project.path.data,'example3'),'001', [project.h,project.w]);
NdepthEst = repmat(imnormalize(depthEst),[1,1,3,1]); %Normalize/add channels for visualization
figure('Name','Input images and estimated depth'),
imshow([imgs(:,:,:,9), imgs(:,:,:,18), imgs(:,:,:,27); ...
        NdepthEst(:,:,:,9), NdepthEst(:,:,:,18), NdepthEst(:,:,:,27)]);
implay([imgs, NdepthEst]); %Play as video

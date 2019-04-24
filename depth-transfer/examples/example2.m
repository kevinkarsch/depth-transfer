%EXAMPLE2 This example shows how to use depth transfer on an input sequence
%         with a rotating camera
%
%%%%%%%%%%%   Begin example2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% Create test data for example2
example2Files = dir( fullfile('example2_data', '*.jpg') ); %Get listing of all images
imgExample2Cell = cell(numel(example2Files),1); %Initialize image cell
for i=1:numel(example2Files) %Read all images as doubles
    imgExample2Cell{i} = im2double(imread(fullfile('example2_data',example2Files(i).name)));
end
imgExample2 = cell2mat(reshape(imgExample2Cell,1,1,1,[])); %Convert cell -> 4D array
%Add data (unless it already exists from a previous run)
if( ~exist(fullfile(project.path.data,'example2','001'), 'dir') )
    createData(fullfile(project.path.data,'example2'), imgExample2); %Reformat data
end

%% Depth transfer
%Set which data to train with
trainFiles = dir(fullfile(project.path.data, 'Make3D-Train*'));
%Infer depth for created data (data/example2/001)
testFile = fullfile('example2', '001');
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
%motionFunc = @segmentGlobalMotionSimple; %We could use this one since it 
                                          % works best on short sequences, 
motionFunc = [];                          % but since there are no moving 
                                          % objects, won't use motion 
                                          % segmentation
%Run depth transfer
depthEst = depthTransfer(project, testFile, trainFiles, depthPrior, motionFunc);

%% Display results
imgs = loadData(fullfile(project.path.data,'example2'),'001', [project.h,project.w]);
NdepthEst = repmat(imnormalize(depthEst),[1,1,3,1]); %Normalize/add channels for visualization
figure('Name','Input images and estimated depth'),
imshow([imgs(:,:,:,1), imgs(:,:,:,4), imgs(:,:,:,7); ...
        NdepthEst(:,:,:,1), NdepthEst(:,:,:,4), NdepthEst(:,:,:,7)]);
implay([imgs, NdepthEst]); %Play as video

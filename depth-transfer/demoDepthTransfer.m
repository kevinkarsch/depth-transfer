%DEMODEPTHTRANSFER This script shows the proper usage of DepthTransfer, 
% and how to create testing and training data
%
EXAMPLES_DIR = 'examples'; %Example directory in root of DepthTransfer
%
%%%%%%%%%%%   Begin demoDepthTransfer   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize a new project (sets up paths, parameters, etc)
h = 460; w = 345; %Inferred depth resolution (output)
Cv = 7; %Number of candidate videos to use for training
Cf = 1; %Number of candidate frames from each video
project = initializeProject(Cv, Cf, [h,w]);

%% Create training data
fprintf('Preparing training data...'); prepTrainTime = tic;
trainFilesDir = fullfile(EXAMPLES_DIR, 'sample_training_data');
trainFiles = dir(fullfile(trainFilesDir, 'img-*.jpg'));
parfor i=1:numel(trainFiles)
    tmpProject = project; %Avoid parfor slicing issue
    [~, name, ~] = fileparts(fullfile(trainFilesDir, trainFiles(i).name));
    basename = name(5:end); %Remove 'img-' prefix
    dataDirName = fullfile(tmpProject.path.data,['Make3D-Train-' basename]);
    if( exist(fullfile(dataDirName, '001'), 'dir') )
        continue; %Training data already exists
    end
    img = imread(fullfile(trainFilesDir, trainFiles(i).name));
    depthFile = dir( fullfile(trainFilesDir, ['depth_sph_corr-' basename '.mat']) );
    foo = load( fullfile(trainFilesDir, depthFile(1).name) );
    depth = foo.Position3DGrid(:,:,4); %Load only depth from laser data
    createData(dataDirName, img, depth, [], false); %false => verbose off
end
fprintf('done. [%6.02fs]\n', toc(prepTrainTime));

%% Create test data for example 1 image (unless it already exists)
img = im2double(imread(fullfile(EXAMPLES_DIR,'demo_data','img-op57-p-016t000.jpg')));
if( ~exist(fullfile(project.path.data,'demo','001'), 'dir') )
    createData(fullfile(project.path.data,'demo'), img);
end

%% Depth transfer
%Set which data to train with
trainFiles = dir(fullfile(project.path.data, 'Make3D-Train*'));
testFile = fullfile('demo', '001');
%Compute prior (average training depth). To save time for future runs, here
% we precompute and store the prior (training data must stay constant).
if( exist(fullfile(EXAMPLES_DIR,'sample_training_prior.mat'), 'file') )
    load(fullfile(EXAMPLES_DIR,'sample_training_prior.mat'));
else
    fprintf('Computing depth prior...'); testTime = tic;
    depthPrior = computePrior(project, trainFiles);
    save(fullfile(EXAMPLES_DIR,'sample_training_prior.mat'), 'depthPrior');
    fprintf('done.   [%6.02fs]\n', toc(testTime));
end
%Set the motion segmentation function here. Since this is a single image, 
% we actially don't need one. See 'examples' directory for using these
% functions, or depthTransfer.m for documentation and available
% segmentation functions
motionFunc = [];
%Run depth transfer
depthEst = depthTransfer(project, testFile, trainFiles, depthPrior, motionFunc);

%% Display results
img = imresize(img,[project.h,project.w]);
NdepthEst = repmat(imnormalize(depthEst),[1,1,3,1]); %Normalize/add channels for visualization
imshow([img, NdepthEst]);

%EVALUATEMAKE3DDATASET Download and evaluate DepthTransfer on the Make3D Depth+Laser dataset
% 
% This script reproduces many of the Make3D results found in 
%  [Karsch et al. ECCV'12], and performs the following tasks:
%    - download the entire Make3D benchmark dataset
%    - reformat the data according to DepthTransfer specifications
%    - estimate depth for each training image in the dataset
%    - evaluate the depth error in several common metrics
% All of these operations can take a *very* long time, so be prepared for a
%  long run. If you have access to Parallel Computing Toolbox, you might
%  want to enable extra threads (i.e. 'matlabpool NUM_PROCESSORS'). With 
%  only one processor on a slower machine, this could take 20 hours!
% The unpacked Make3D dataset is almost 1GB, and the script produces 
%  another 3GB of data, so please make sure >4GB of space is available on 
%  your disk.
%
% NOTE: If you already have the Make3D dataset in its original format, 
%  modify the below directory pointers to their correct paths (relative to 
%  this directory or absolute). Otherwise, these don't need to be modified.
MAKE3D_TRAIN_IMG_DIR = 'Train400Img';
MAKE3D_TRAIN_DEPTH_DIR = 'Train400Depth';
MAKE3D_TEST_IMG_DIR = 'Test134';
MAKE3D_TEST_DEPTH_DIR = 'Gridlaserdata';
%
%%%%%%%%%%%   Begin evaluateMake3dDataset   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Download and/or initialize Make3D Laser+Image data
% This script requires the Make3D range image data, found here: 
%   http://make3d.cs.cornell.edu/data.html
%
% Use the code below to download the dataset and unpack in the current
% directory. If you already have the dataset or prefer to download/unpack
% yourself, comment the following lines (or set DOWNLOAD_MAKE3D=false), and
% ensure that the directory pointers (MAKE3D_*_DIR) are set properly above.
if( ~exist(MAKE3D_TRAIN_IMG_DIR, 'dir') )
    fprintf('Downloading Make3D training images...');
    untar('http://cs.stanford.edu/group/reconstruction3d/Train400Img.tar.gz', MAKE3D_TRAIN_IMG_DIR);
    fprintf('done.\n');
end
if( ~exist(MAKE3D_TRAIN_DEPTH_DIR, 'dir') )
    fprintf('Downloading Make3D training depth maps...');
    untar('http://cs.stanford.edu/group/reconstruction3d/Train400Depth.tgz', fileparts(MAKE3D_TRAIN_DEPTH_DIR));
    fprintf('done.\n');
end
if( ~exist(MAKE3D_TEST_IMG_DIR, 'dir') )
    fprintf('Downloading Make3D testing images...');
    untar('http://www.cs.cornell.edu/~asaxena/learningdepth/Test134.tar.gz', fileparts(MAKE3D_TEST_IMG_DIR));
    fprintf('done.\n');
end
if( ~exist(MAKE3D_TEST_DEPTH_DIR, 'dir') )
    fprintf('Downloading Make3D testing depth maps...');
    untar('http://www.cs.cornell.edu/~asaxena/learningdepth/Test134Depth.tar.gz', fileparts(MAKE3D_TEST_DEPTH_DIR));
    fprintf('done.\n');
end

%% Set project parameters
%Set project/output size. All images and depths will be scaled to match 
% these dimensions. (These were chosen to preserve the aspect ratio of the
% Make3D color images without significantly stretching the depth maps 
% horizontally).
h = 460; w = 345;
Cv = 7; %Number of candidate videos to use for training
Cf = 1; %Number of candidate frames from each video (since this code 
        % demonstrates on a single-image dataset, each video only contains
        % one frame).
project = initializeProject(Cv, Cf, [h,w]);

%% Reformat Make3D data into folders/structs that DepthTransfer code needs
% Note: once the Make3D data has been reformatted, the original data
%   (MAKE3D_*_DIR directories) can be removed.
% Add training data
fprintf('Preparing training data...'); prepTrainTime = tic;
trainM3DImgFiles = dir(fullfile(MAKE3D_TRAIN_IMG_DIR, 'img-*.jpg'));
trainM3DDepthFiles = dir(fullfile(MAKE3D_TRAIN_DEPTH_DIR, 'depth_sph_corr-*.mat'));
assert(numel(trainM3DImgFiles)==400, ...
    'Incorrect number of Make3D training images (usually due to download/unpacking error)\n');
assert(numel(trainM3DDepthFiles)==400, ...
    'Incorrect number of Make3D training depths (usually due to download/unpacking error)\n');
parfor i=1:numel(trainM3DImgFiles)
    tmpProject = project; %Avoid parfor slicing issue
    [~, name, ~] = fileparts(fullfile(MAKE3D_TRAIN_IMG_DIR, trainM3DImgFiles(i).name));
    basename = name(5:end); %Remove 'img-' prefix
    dataDirName = fullfile(tmpProject.path.data,['Make3D-Train-' basename]);
    if( exist(fullfile(dataDirName, '001', 'features.mat'), 'file') )
        continue; %Data already processed
    end
    img = imread(fullfile(MAKE3D_TRAIN_IMG_DIR, trainM3DImgFiles(i).name));
    depthFile = dir( fullfile(MAKE3D_TRAIN_DEPTH_DIR, ['depth_sph_corr-' basename '.mat']) );
    foo = load( fullfile(MAKE3D_TRAIN_DEPTH_DIR,depthFile(1).name) );
    depth = foo.Position3DGrid(:,:,4); %Load only depth from laser data
    createData(dataDirName, img, depth, [], false); %false => verbose off
end
fprintf('done. [%6.02fs]\n', toc(prepTrainTime));

% Add testing data
fprintf('Preparing testing data...'); prepTestTime = tic;
testM3DImgFiles = dir(fullfile(MAKE3D_TEST_IMG_DIR, 'img-*.jpg'));
testM3DDepthFiles = dir(fullfile(MAKE3D_TEST_DEPTH_DIR, 'depth_sph_corr-*.mat'));
assert(numel(testM3DImgFiles)==134, ...
    'Incorrect number of Make3D test images (usually due to download/unpacking error)\n');
assert(numel(testM3DDepthFiles)==134, ...
    'Incorrect number of Make3D test depths (usually due to download/unpacking error)\n');
testImg = zeros([project.h, project.w, 3, 134]);
parfor i=1:numel(testM3DImgFiles)
    tmpProject = project; %Avoid parfor slicing issue
    [~, name, ~] = fileparts(fullfile(MAKE3D_TEST_IMG_DIR, testM3DImgFiles(i).name));
    basename = name(5:end); %Remove 'img-' prefix
    dataDirName = fullfile(tmpProject.path.data,['Make3D-Test-' name]);
    img = imread(fullfile(MAKE3D_TEST_IMG_DIR, testM3DImgFiles(i).name));
    testImg(:,:,:,i) = im2double(imresize(img, [project.h, project.w]));
    if( exist(fullfile(dataDirName, '001', 'features.mat'), 'file') )
        continue; %Data already processed
    end
    depthFile = dir( fullfile(MAKE3D_TEST_DEPTH_DIR, ['depth_sph_corr-' basename '.mat']) );
    foo = load( fullfile(MAKE3D_TEST_DEPTH_DIR,depthFile(1).name) );
    depth = foo.Position3DGrid(:,:,4); %Load only depth from laser data
    createData(dataDirName, img, depth, [], false); %false => verbose off
end
fprintf('done.  [%6.02fs]\n', toc(prepTestTime));

%% Run depth transfer on each of the Make3D testing images
%Set train/test split
trainFiles = dir(fullfile(project.path.data, 'Make3D-Train*'));
testFiles = dir(fullfile(project.path.data, 'Make3D-Test*'));

%Compute prior (average training depth)
fprintf('Computing depth prior...'); testTime = tic;
depthPrior = computePrior(project, trainFiles);
fprintf('done.   [%6.02fs]\n', toc(testTime));

%Run depth transfer
fprintf(['\nRunning depth transfer on all test images\n', ...
           '=========================================\n']); testTime = tic;
depthEst = zeros(h,w,numel(testFiles)); %Estimated depth from DepthTransfer
depthTrue = zeros(h,w,numel(testFiles)); %Ground truth depth from dataset
for i=1:numel(testFiles)
    fprintf('\nProcessing %d of %d (%s)\n\n', i, numel(testFiles), testFiles(i).name);
    testClip = fullfile(testFiles(i).name, '001');
    [depthEst(:,:,i) depthTrue(:,:,i)] = ...
        depthTransfer(project, testClip, trainFiles, depthPrior, []);
    save('RESULT.mat', 'depthEst', 'depthTrue');
end
fprintf('done.   [%6.02fs]\n', toc(testTime));

%% Display and evaluate results
allDepthsLog10 = imnormalize(permute([depthEst, depthTrue],[1,2,4,3]));
implay([testImg, repmat(allDepthsLog10,[1,1,3,1])]);

rel = reshape( mean(mean( abs(depthEst-depthTrue)./depthTrue, 2),1), [],1);
lg10 = reshape( mean(mean( abs(log10(depthEst)-log10(depthTrue)), 2),1), [],1);
rmse = reshape( sqrt(mean(mean( abs(depthEst-depthTrue).^2, 2),1)), [],1);
fprintf(['\nError averaged over all test data\n', ...
           '=================================\n', ...
           '%10s %10s %10s\n%10.4f %10.4f %10.4f\n'], ...
           'relative', 'log10', 'rmse', mean(rel), mean(lg10), mean(rmse));
		   
%If everything ran properly, the final line of output should say:
% Error averaged over all test data
% =================================
%   relative      log10       rmse
%     0.3608     0.1480    15.1446
%
%You will also see a figure displaying each frame's result (from left to 
% right: image, estimated depth, ground truth depth). Depth images are in 
% log10 space and normalized over the entire dataset for visualization.

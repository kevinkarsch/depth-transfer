function depth = depthTransferImg(img, trainClips, prior, project)
%DEPTHTRANSFER Estimate depth for an image/video using the DepthTransfer approach
%
% Input:
%  img        - Input image for which to estimate depth
%  trainClips - Struct array such that trainClips(i).name contains the 
%               i^{th} training data directory (relative to
%               project.path.data). MATLAB's built in dir() function can 
%               create these structs; ex:
%                   trainClips = ...
%                       dir(fullfile(project.path.data, [BASE_NAME '*']);
%               where BASE_NAME is a name shared by all training 
%               directories
%  prior      - Depth prior used during optimization; usually set as the
%               average of all training depths, which can be computed with
%               computePrior(...). Size should be [project.h, project.w].
%               If prior=[] or not provided, then no prior is used.
%
% Output:
%  depthEst  - Estimated depth
%
%%%%%%%%%%%   Begin depthTransferImg   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
[h,w,~] = size(img);
if(~exist('project', 'var') || isempty(project))
    project = initializeProject(7, 1, [h,w]);
    project.para.depthTransferParams.smoothCoeff_tem = 0;
    project.para.depthTransferParams.motionCoeff     = 0;
end
if(~exist('prior','var') || isempty(prior))
    prior = zeros(h,w);
    project.para.depthTransferParams.priorCoeff      = 0;
end
K = 1;
imgs = im2double(imresize(img,[h,w],'bilinear'));
flows = zeros(h,w,2);
features = computeFeatures(img, flows);
videogist = features.videogist;
flowgist = features.flowgist;
gist = features.gist;
motionseg = zeros(h,w,1);

%% Load training features (or compute them if they don't exist)
% This can take some time on the first run, but features will only be
% computed ONCE for each testing data and stored on disk. Any time these 
% data are used again, the features will be loaded loaded from disk
nntime = tic; fprintf('Precomputing features for candidate selection...');
training.N = numel(trainClips);
training.videogist = cell(training.N, 1);
training.flowgist = cell(training.N, 1);
training.gist = cell(training.N, 1);
for i=1:training.N
    trainDataDir = fullfile(project.path.data, trainClips(i).name);
    foo = load(fullfile(trainDataDir, 'info.mat'));
    trainClipInfo = foo.clips;
    for j=1:numel(trainClipInfo)
        trainClipName = fullfile(trainDataDir, trainClipInfo{j}.name);
        if( ~exist(fullfile(trainClipName, 'features.mat'), 'file') )
            [ftimgs, ~, ftflows] = loadData(trainDataDir, trainClipInfo{j}.name);            
            features = computeFeatures(ftimgs, ftflows); %#ok<NASGU>
            save(fullfile(trainClipName, 'features.mat'), 'features');
        end
        foo = load(fullfile(trainClipName, 'features.mat'));
        training.videogist{i}{j} = foo.features.videogist;
        training.flowgist{i}{j} = foo.features.flowgist;
        training.gist{i}{j} = foo.features.gist;
    end
end
fprintf('done. [%6.2fs]\n', toc(nntime));

%% Find Nearest neighbors and SIFT flow warping operator
sftime = tic; fprintf('Computing candidates and siftflow operators...\n');
%Compute candidate videos
candidates = cell(K,1);
videogistError = zeros(training.N,1);
for i=1:training.N
    videogistError(i) = min( cellfun(@(x) sum((x-videogist).^2), training.videogist{i}) );
end
[~,idx] = sort(videogistError, 'ascend');
train_videoidx = idx(1:project.Cv);

%Compute candidate frames (project.Cf frames selected from candidate videos)
ftsWeight = project.para.depthTransferParams.gf_mix;
for frame=1:K
    tmpProject = project; %Avoid parfor slicing issue
    tmpTraining = training; %Avoid parfor slicing issue
    tmpTrainClips = trainClips; %Avoid parfor slicing issue
    tmpTrain_videoidx = train_videoidx; %Avoid parfor slicing issue
    for i=1:tmpProject.Cv
        j = tmpTrain_videoidx(i);
        trainDataDir = fullfile(tmpProject.path.data, tmpTrainClips(j).name);
        foo = load(fullfile(trainDataDir, 'info.mat'));
        trainClipInfoj = foo.clips;
        numTrainClipsj = numel(trainClipInfoj);
        allerrsj = cell(numTrainClipsj,1);
        clipidx_listj = cell(numTrainClipsj,1);
        frameidx_listj = cell(numTrainClipsj,1);
        for ii=1:numTrainClipsj
            trflowgist = tmpTraining.flowgist{j}{ii};
            flowgisterr = repmat(flowgist(frame,:),[size(trflowgist,1),1])-trflowgist;
            trgist = tmpTraining.gist{j}{ii};
            gisterr = repmat(gist(frame,:),[size(trgist,1),1])-trgist;
            allerrsj{ii} = (1-ftsWeight).*sum(flowgisterr.^2,2) + ftsWeight.*sum(gisterr.^2,2);
            clipidx_listj{ii} = ii.*ones(numel(allerrsj{ii}),1); %Which clip to candidate is from
            frameidx_listj{ii} = (1:numel(allerrsj{ii}))'; %Which frame of the given clip
        end
        [~,idx] = sort(cell2mat(allerrsj),'ascend'); %Sort all frames by error
        clipidx_listj = cell2mat(clipidx_listj); %Vectorize
        frameidx_listj = cell2mat(frameidx_listj); %Vectorize
        candidates{frame}{i}.videodir = tmpTrainClips(j).name; %Get training data name
        for k=1:tmpProject.Cf
            candidates{frame}{i}.subvideodirs{k} = trainClipInfoj{clipidx_listj(idx(k))}.name; %Get clip name
            candidates{frame}{i}.clipInfo{k} = trainClipInfoj{clipidx_listj(idx(k))};          %and info
        end
        candidates{frame}{i}.frameidx = frameidx_listj(idx(1:tmpProject.Cf));
    end
end

% Find warping operator with SIFTflow
sfwarpx = zeros(h,w,project.Cv*project.Cf,K);
sfwarpy = zeros(h,w,project.Cv*project.Cf,K);
SIFTError = zeros(h,w,project.Cv*project.Cf,K);
for frame=1:K
    sift = mexDenseSIFT(imgs(:,:,:,frame), project.para.SIFT.cellsize, project.para.SIFT.gridspacing);
    sfwarpxi = zeros(h,w,project.Cv*project.Cf);
    sfwarpyi = zeros(h,w,project.Cv*project.Cf);
    SIFTErrori = zeros(h,w,project.Cv*project.Cf);
    fprintf('\tFrame %d:\n', frame);
    parfor idx=1:project.Cv*project.Cf
        tmpProject = project; %Avoid parfor slicing issue
        tmpCandidates = candidates; %Avoid parfor slicing issue
        i = floor((idx-1)./tmpProject.Cf)+1;
        j = mod(idx-1,tmpProject.Cf)+1;
        train_dir = fullfile(project.path.data, tmpCandidates{frame}{i}.videodir);
        train_subdir = fullfile(train_dir, tmpCandidates{frame}{i}.subvideodirs{j});
        train_imgdir = tmpCandidates{frame}{i}.clipInfo{j}.imgDir;
        train_frames = dir(fullfile(train_subdir, train_imgdir,'*.png'));
        train_frames = train_frames(tmpCandidates{frame}{i}.frameidx);
        train_imgfile = fullfile(train_subdir, train_imgdir, train_frames(j).name);
        train_img = imresize( im2double(imread( train_imgfile )), [h,w]);
        train_sift = mexDenseSIFT(train_img, project.para.SIFT.cellsize, project.para.SIFT.gridspacing);  
        fprintf('\t\tCandidate %d: %s\n', idx, fullfile(train_subdir, train_imgdir, train_frames(j).name));
        [sfwarpxi(:,:,idx), sfwarpyi(:,:,idx), ~] = ...
            SIFTflowc2f(sift, train_sift, project.para.SIFTflowpara);
        SIFTErrori(:,:,idx) = sum(abs( double(sift)./255 - ...
            double(mexWarpImageInt(sift,train_sift,sfwarpxi(:,:,idx), sfwarpyi(:,:,idx)))), 3);
    end
    sfwarpx(:,:,:,frame) = sfwarpxi;
    sfwarpy(:,:,:,frame) = sfwarpyi;
    SIFTError(:,:,:,frame) = SIFTErrori;
end
fprintf('done. [%6.2fs]\n', toc(sftime));

%% Warp candidate data to input video
candidatewarptime = tic; fprintf('Warping candidates to test image(s)...');
Depth       = zeros(h,w,project.Cv*project.Cf,K);
DepthX      = zeros(h,w,project.Cv*project.Cf,K);
DepthY      = zeros(h,w,project.Cv*project.Cf,K);
FlowError   = zeros(h,w,project.Cv*project.Cf,K);
Mask        = zeros(h,w,project.Cv*project.Cf,K);
for frame=1:K
    flowmag      = sqrt(sum(flows(:,:,:,frame).^2,3));
    Depthi       = zeros(h,w,project.Cv*project.Cf);
    DepthXi      = zeros(h,w,project.Cv*project.Cf);
    DepthYi      = zeros(h,w,project.Cv*project.Cf);
    FlowErrori   = zeros(h,w,project.Cv*project.Cf);
    Maski        = zeros(h,w,project.Cv*project.Cf);
    parfor idx=1:project.Cv*project.Cf
        tmpProject = project; %Avoid parfor slicing issue
        tmpCandidates = candidates; %Avoid parfor slicing issue
        %Find correct frame
        i = floor((idx-1)./tmpProject.Cf)+1;
        j = mod(idx-1,tmpProject.Cf)+1;
        train_info = tmpCandidates{frame}{i}.clipInfo{j};
        train_dir = fullfile(tmpProject.path.data, tmpCandidates{frame}{i}.videodir);
        train_subdir = fullfile(train_dir, tmpCandidates{frame}{i}.subvideodirs{j});
        train_frames = dir(fullfile(train_subdir, train_info.imgDir, '*.png'));
        train_frames = train_frames(tmpCandidates{frame}{i}.frameidx);
        %Load depth, flow, and hole masks. Resize to project dimensions
        train_depth = im2double(imread(fullfile(train_subdir, train_info.depthDir, train_frames(j).name)));
        train_depth = train_depth.*(train_info.depth_bounds(2)-train_info.depth_bounds(1)) + train_info.depth_bounds(1);
        train_depth = imresize(train_depth, [h,w]);
        train_depthx = imfilter(train_depth, [-1,1,0]);
        train_depthy = imfilter(train_depth, [-1,1,0]');
        train_flow = im2double(imread(fullfile(train_subdir, train_info.flowDir, train_frames(j).name)));
        train_flow = train_flow(:,:,1:2,:).*(train_info.flow_bounds(2)-train_info.flow_bounds(1)) + train_info.flow_bounds(1);
        flow_rescale = permute([h,w]./[size(train_flow,1),size(train_flow,2)],[1,3,2]);
        train_flow = imresize(train_flow, [h,w], 'bilinear').*repmat(flow_rescale,[h,w,1]);
        train_mask = im2double(imread(fullfile(train_subdir, train_info.maskDir, train_frames(j).name)));
        train_mask = imresize(train_mask, [h,w], 'bilinear');
        %warp and store
        vx = sfwarpx(:,:,idx,frame);
        vy = sfwarpy(:,:,idx,frame);
        Depthi(:,:,idx) = mexWarpImageInt(train_depth, train_depth, vx, vy);
        DepthXi(:,:,idx) = mexWarpImageInt(train_depthx, train_depthx, vx, vy);
        DepthYi(:,:,idx) = mexWarpImageInt(train_depthy, train_depthy, vx, vy);
        FlowErrori(:,:,idx) = abs(double(mexWarpImageInt(flowmag,sqrt(sum(train_flow.^2,3)),vx,vy)) - flowmag);
        Maski(:,:,idx) = mexWarpImageInt(train_mask,train_mask, vx, vy);
    end
    Depth(:,:,:,frame) = Depthi;
    DepthX(:,:,:,frame) = DepthXi;
    DepthY(:,:,:,frame) = DepthYi;
    FlowError(:,:,:,frame) = FlowErrori;
    Mask(:,:,:,frame) = Maski;
end
fprintf('done. [%6.2fs]\n', toc(candidatewarptime));

%% Infer depth with label transfer
fprintf('Inferring depth...\n');
%Setup prior
if( ~exist('prior','var') || isempty(prior) )
    prior = zeros(h,w);
    project.para.depthTransferParams.priorCoeff = 0;
end
Prior = repmat(imresize(prior, [h,w], 'bilinear'),[1,1,1,K]);
%Run inference
depth = depthVideoInference(imgs, flows, motionseg, Depth, DepthX, DepthY, SIFTError, FlowError, Mask, Prior, project.para.depthTransferParams);

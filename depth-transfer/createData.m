function createData( dataDirectory, img, depth, holeMask, verbose )
%CREATEDATA Reformats and saves data that DepthTransfer code expects
%
% Input:
%  dataDirectory - directory to save data (img, depth, etc) to [string]
%  img           - 4D img/video data [height x width x {1,3} x numFrames] 
%                  OR a cell array of 4D video clips (use cell format if 
%                  video should be split into smaller clips)
%  depth(=[])    - 4D depth data [height x width x 1 x numFrames] OR a
%                  cell array of 4D depth clips. Can also be formatted as 
%                  3D (without singleton dim). img and depth do NOT need to
%                  have the same height or width (but numFrames should be 
%                  the same though)
%  holeMask(=[]) - 4D binary hole data [height x width x 1 x numFrames] OR 
%                  a cell array of 4D hole info (should be consistent with
%                  depth). Use this variable if depth contains known pixels 
%                  of missing information (pixels == 0 will be 
%                  interpolated), otherwise disregard it or leave it as 
%                  empty ([]). Can also be formatted as 3D (without 
%                  singleton dim)
%  verbose(=true)- Print timing/debug information (default is true)
%    
% Note: To achieve temporal consistency and motion estimation, optical 
%   flow must be computed between neighboring frames. We use the publicly 
%   available optical flow code from Ce Liu, found here:
%   http://people.csail.mit.edu/celiu/OpticalFlow/
%   opticalflow.m (found in this directory) provides an interface for Liu's
%   code. To use any other optical flow module, either edit opticalflow.m,
%   or modify opticalFlowFunc to point to your own flow function (make sure
%   the output data format is the same as in opticalflow.m).
%   Set opticalFlowFunc = [] if you do not need flow (i.e. running on
%   single images).
opticalFlowFunc = @opticalflow;
%
%%%%%%%%%%%   Begin createData   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Convert everything to cell arrays
    if( ~iscell(img) )
        img = {img};
        if( ~exist('depth', 'var') )
            depth = {[]};
        else
            depth = {depth};
        end
        if( ~exist('holeMask', 'var') )
            holeMask = {[]};
        else
            holeMask = {holeMask};
        end
    end
    if( ~exist('depth', 'var') )
        depth = repmat({[]}, size(img));
    end
    if( ~exist('holeMask', 'var') )
        holeMask = repmat({[]}, size(img));
    end
    if( ~exist('verbose', 'var') )
        verbose = true;
    end
    
    %Input validation
    [h,w,d,~] = size(img{1});
    [h2,w2,~,~] = size(depth{1});
    for i=1:numel(img)
        [hi,wi,~,Ki] = size(img{i});
        assert(h==hi && w==wi, 'Input dimension mismatch (clip #%03d)\n', i);
        if( isempty(depth{i}) && isempty(holeMask{i})  )
            depth{i} = zeros(h,w,Ki);
            holeMask{i} = true(h,w,Ki);
        elseif( isempty(holeMask{i}) )
            holeMask{i} = true(h2,w2,Ki);
        else
            try
                depth{i} = reshape(depth{i}, [h2,w2,Ki]);
                holeMask{i} = reshape(holeMask{i}, [h2,w2,Ki]);
            catch %#ok<CTCH>
                error('depth and holeMask must have same dimensions (clip #%03d)\n', i);
            end
        end
    end
    
    %Create output dir
    if( ~exist(dataDirectory,'dir') )
        mkdir(dataDirectory);
        clips = [];
        save(fullfile(dataDirectory,'info.mat'), 'clips');
    else
        warning('%s already exists. Adding data as additional clip(s).\n', dataDirectory); %#ok<WNTAG>
        foo = load(fullfile(dataDirectory,'info.mat'));
        clips = foo.clips;
    end
    
    %Save data
    if(verbose), fprintf('Creating DepthTransfer data at: %s\n', dataDirectory); end
    for i=1:numel(img)
        nc = numel(clips) + 1;
        if(verbose), fprintf('Processing clip %03d\n', nc); end
        K = size(img{i},4);
        
        %Compute optical flow (if necessary)
        if(verbose), fprintf('\tComputing optical flow...'); flowtime = tic; end
        flow = zeros(h,w,2,K);
        if( ~isempty(opticalFlowFunc) )
            tmpimg = img{i};
            tmpimg_next = img{i}(:,:,:,2:K);
            parfor j=1:K-1
                flow(:,:,:,j) = opticalFlowFunc(tmpimg(:,:,:,j), tmpimg_next(:,:,:,j));
            end
        end
        if(verbose), fprintf('done. [%6.02fs]\n', toc(flowtime)); end
        
        %Fill depth holes (if necessary)
        if( any(holeMask{i}(:)==0) )
            if(verbose), fprintf('\tFilling depth holes...'); filltime = tic; end
            img_resize = reshape(imresize(img{i}, [h2,w2], 'bilinear'), [h2,w2,d,K]);
            depth_filled = fillDepthHoles(img_resize, depth{i}, flow, holeMask{i});
            if(verbose), fprintf('done. [%6.02fs]\n', toc(filltime)); end
        else
            depth_filled = depth{i};
        end
        
        %Write data
        if(verbose), fprintf('\tSaving data...'); savetime = tic; end
        clipInfo.name = sprintf('%03d', nc);
        clipInfo.size = [h,w];
        clipInfo.numFrames = K;
        clipInfo.imgDir = 'img';
        clipInfo.depth0Dir = 'depth0';
        clipInfo.depthDir = 'depth';
        clipInfo.flowDir = 'flow';
        clipInfo.maskDir = 'mask';
        clipInfo.flow_bounds = [min(flow(:)), max(flow(:))];
        clipInfo.depth0_bounds = [min(depth{i}(:)), max(depth{i}(:))];
        clipInfo.depth_bounds = [min(depth_filled(:)), max(depth_filled(:))];
        clips{nc} = clipInfo; %#ok<AGROW>
        mkdir(fullfile(dataDirectory, clipInfo.name));
        mkdir(fullfile(dataDirectory, clipInfo.name, clipInfo.imgDir));
        mkdir(fullfile(dataDirectory, clipInfo.name, clipInfo.depth0Dir));
        mkdir(fullfile(dataDirectory, clipInfo.name, clipInfo.depthDir));
        mkdir(fullfile(dataDirectory, clipInfo.name, clipInfo.flowDir));
        mkdir(fullfile(dataDirectory, clipInfo.name, clipInfo.maskDir));
        %Convert floating point images to 16 bit uints
        depth0 = uint16(round(65535.*imnormalize(depth{i})));
        depth_filled = uint16(round(65535.*imnormalize(depth_filled)));
        flow3 = cat(3, imnormalize(flow), zeros(h,w,1,K));
        flow3 = uint16(round(65535.*flow3));
        %Make sure holeMask is stored as 1 bit
        holeMask{i} = logical(holeMask{i}~=0);
        for j=1:K
            namej = sprintf('%04d.png', j-1);
            imwrite(img{i}(:,:,:,j), fullfile(dataDirectory, clipInfo.name, clipInfo.imgDir, namej));
            imwrite(depth0(:,:,j), fullfile(dataDirectory, clipInfo.name, clipInfo.depth0Dir, namej));
            imwrite(depth_filled(:,:,j), fullfile(dataDirectory, clipInfo.name, clipInfo.depthDir, namej));
            imwrite(flow3(:,:,:,j), fullfile(dataDirectory, clipInfo.name, clipInfo.flowDir, namej));
            imwrite(holeMask{i}(:,:,j), fullfile(dataDirectory, clipInfo.name, clipInfo.maskDir, namej));
        end
        %Save info struct
        save(fullfile(dataDirectory,'info.mat'), 'clips');
        if(verbose), fprintf('done. [%6.02fs]\n', toc(savetime)); end
    end
end


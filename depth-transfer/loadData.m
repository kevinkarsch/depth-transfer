function [img, depth, flow, holeMask, depth0, clipInfo, features] = ...
    loadData( dataDirectory, clipIndicesOrName, resizeDims )
%LOADDATA Loads a directory containing video clips of DepthTransfer data
%
% Input:
%  dataDirectory          - Directory containing DepthTransfer formatted 
%                           data
%  clipIndicesOrName(=[]) - Either the clip index(es) or the name of the 
%                           clip within dataDirectory that should be 
%                           loaded
%  resizeDims(=[])        - 1x2 integer spatial dimensions that all data 
%                           (img, depth, flow, etc) should be rescaled to
%
% Output:
%  img      - 4D array [height x width x {1,3} x numFrames] of image/video 
%             data. If clipIndicesOrName is a vector, then img will be a
%             cell of 4D arrays of size numel(clipIndicesOrName)
%  depth    - Loaded depth data. Of size [height x width x 1 x numFrames]
%  flow     - Loaded optical flow data. Of size 
%             [height x width x 2 x numFrames]
%  holeMask - Loaded hole mask data. holeMask(i)==0 implies depth pixel i
%             is unknown or interpolated. Same size as depth
%  depth0   - Loaded original depth data. The difference between depth and
%             depth0 is depth has had its holes interpolated, depth0 has 
%             not. If there are no holes (ie all(holeMask(:))==1), then
%             depth and depth0 are the same. Same size as depth
%  clipInfo - Info struct about the loaded clip (or cell of info structs if
%             clipIndicesOrName is a vector)
%  features - Precomputed features used for selecting candidate
%             video/frames in depthTransfer(...). If features have not been
%             precomputed, features=[].
%
%%%%%%%%%%%   Begin loadData   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    foo = load(fullfile(dataDirectory,'info.mat'));
    clipInfo = foo.clips;
    if( ~exist('clipIndicesOrName','var') || isempty(clipIndicesOrName) )
        clipIndex = 1:numel(clipInfo);
    elseif( ischar(clipIndicesOrName) )
        clipIndex = find(cellfun(@(x) strcmp(x.name,clipIndicesOrName), clipInfo));
    else
        clipIndex = clipIndicesOrName;
    end
    loadFeatures = (nargout>6);
    
    %Load all clips
    img = cell(numel(clipIndex),1);
    holeMask = cell(numel(clipIndex),1);
    depth = cell(numel(clipIndex),1);
    depth0 = cell(numel(clipIndex),1);
    flow = cell(numel(clipIndex),1);
    features = cell(numel(clipIndex),1);
    i = 1;
    for nc=clipIndex
        datapath = fullfile(dataDirectory, clipInfo{nc}.name);
        %Load images
        img{i} = loadPNGDir(fullfile(datapath, clipInfo{nc}.imgDir), clipInfo{nc});
        holeMask{i} = loadPNGDir(fullfile(datapath, clipInfo{nc}.maskDir), clipInfo{nc});
        depth{i} = loadPNGDir(fullfile(datapath, clipInfo{nc}.depthDir), clipInfo{nc});
        depth{i} = depth{i}.*(clipInfo{nc}.depth_bounds(2)-clipInfo{nc}.depth_bounds(1)) + ...
            clipInfo{nc}.depth_bounds(1); %Rescale
        depth0{i} = loadPNGDir(fullfile(datapath, clipInfo{nc}.depth0Dir), clipInfo{nc});
        depth0{i} = depth0{i}.*(clipInfo{nc}.depth0_bounds(2)-clipInfo{nc}.depth0_bounds(1)) + ...
            clipInfo{nc}.depth0_bounds(1); %Rescale
        flow{i} = loadPNGDir(fullfile(datapath, clipInfo{nc}.flowDir), clipInfo{nc});
        flow{i} = flow{i}.*(clipInfo{nc}.flow_bounds(2)-clipInfo{nc}.flow_bounds(1)) + ...
            clipInfo{nc}.flow_bounds(1); %Rescale
        %Resize if necessary
        if( exist('resizeDims','var') && ~isempty(resizeDims) )
            img{i} = imresize3(img{i}, resizeDims);
            holeMask{i} = double(imresize3(holeMask{i}, resizeDims)>0.5);
            depth{i} = imresize3(depth{i}, resizeDims);
            depth0{i} = imresize3(depth0{i}, resizeDims);
            flow{i} = imresize3(flow{i}, resizeDims);
            %Rescale flow vectors based on new size, remove empty dimension
            scale = resizeDims./clipInfo{nc}.size;
            flow{i} = cat(3, flow{i}(:,:,1,:).*scale(1),  flow{i}(:,:,2,:).*scale(2));
            clipInfo{nc}.size = resizeDims;
        else
            flow{i}(:,:,3,:) = []; %Remove empty dimension from flow
        end
        %Load features if output argument is provided
        featureFile = fullfile(dataDirectory, clipInfo{nc}.name, 'features.mat');
        if(loadFeatures && exist(featureFile, 'file'))
            foo = load(featureFile); 
            features{i} = foo.features;
        else
            features{i} = [];
        end
        i=i+1;
    end
    
    %Remove from cell array if only one cell
    if(numel(clipIndex)==1) 
        img = img{1}; 
        depth = depth{1}; 
        flow = flow{1}; 
        holeMask = holeMask{1}; 
        depth0 = depth0{1};
        clipInfo = clipInfo{1};
        if(loadFeatures)
            features = features{1};
        end
    end
end

%Loads a directory of PNGs given the directory and clipInfo
function pngs = loadPNGDir(imgdir, clipInfo)
    filenames = dir(fullfile(imgdir, '*.png'));
    %Error checking
    if(numel(filenames)<1)
        pngs = [];
        return;
    end
    [h,w,d] = size(imread( fullfile(imgdir,filenames(1).name )));
    assert(numel(filenames) == clipInfo.numFrames, ...
        'Clip %s: numFrames does not match # of pngs in %s\n', clipInfo.name, imgdir);
    %Read images
    pngs = zeros([h, w, d, clipInfo.numFrames]);
    for i=1:clipInfo.numFrames
        pngs(:,:,:,i) = im2double(imread( fullfile(imgdir,filenames(i).name) ));
    end
end

%Generalizes imresize to videos (as well as images)
function rI = imresize3(I, s)
    rI = reshape(imresize(I, s, 'bilinear'), [s, size(I,3), size(I,4)]);
end

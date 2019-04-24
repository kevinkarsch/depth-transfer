function features = computeFeatures(img, flow) 
%COMPUTEFEATURES Computes feature descriptors used for DepthTransfer candidate matching
% Input:
%  img  - 4D video/img data [height x width x {1,3} x numFrames]. Can also
%         be a cell array of 4D data.
%  flow - 4D optical flow data [height x width x 2 x numFrames]. Can also
%         be a cell array of 4D data.
% Output:
%  features - Struct containing descriptors (cell of structs if img is a 
%             cell array). See corresponding feature functions for more 
%             details (video2gist, im2gist, flow2flowgist).
%%%%%%%%%%%   Begin computeFeatures   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(~iscell(img)) %Convert to cell for looping
        img = {img};
        flow = {flow};
    end
    features = cell(numel(img),1);
    for i=1:numel(img)
        [~,~,~,K] = size(img{i});
        %GIST descriptor for entire video
        features{i}.videogist = video2gist(img{i}); 
        %Flow descriptor per frame
        features{i}.flowgist = zeros(K, flow2flowgist());
        for j=1:K
            features{i}.flowgist(j,:) = flow2flowgist(flow{i}(:,:,:,j));
        end
        %GIST descriptor per frame
        features{i}.gist = zeros(K, im2gist());
        for j=1:K
            features{i}.gist(j,:) = im2gist(img{i}(:,:,:,j));
        end
    end
    if(numel(features)==1)
        features = features{1};
    end
end

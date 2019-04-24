function M = segmentGlobalMotionHeuristic(imgs, flows) %#ok<INUSD>
%SEGMENTGLOBALMOTIONHEURISTIC Performs a binary segmentation of moving pixels
% Uses image stitching and global motion estimation to detect moving pixels
% for rotating and/or variable focal length image sequences. The final
% segmentation comes from a combination of flow and image pixel values.
%
% Input:
%  imgs  - Image/video data of size [height x width x {1,3} x numFrames]
%  flows - Optical flow data corresponding to imgs
% Output:
%  M - Binary motion segmentation of estimated moving pixels in imgs of 
%          size [height x width x 1 x numFrames]
%%%%%%%%%%%   Begin segmentGlobalMotionHeuristic   %%%%%%%%%%%%%%%%%%%%%%%%
    K = size(imgs,4);
    
    %Stitch a background mosaic
    samples = min(5, K); %Only use a few of the images, 5 or K, whichever is less
    bg = stitchImages(imgs(:,:,:,round(linspace(1,K,samples))));
    
    %Warp the input images to the background
    [h,w,d] = size(bg);
    warp = zeros(h,w,d,K);
    H = cell(K,1);
    parfor i=1:K
        [warp(:,:,:,i), H{i}] = registerImagesHomography(bg, imgs(:,:,:,i));
    end
    
    %Compute optical flow on the warped image (for weighting)
    warpnext = warp(:,:,:,2:end);
    flows = zeros(h,w,2,K);
    parfor i=1:K-1
        flows(:,:,:,i) = opticalflow(warp(:,:,:,i), warpnext(:,:,:,i));
    end

    %Motion segmentation heuristic (same as in segmentMotionHeuristic.m)
    moseg = sqrt(sum(flows.^2,3)).*mean((warp-repmat(bg,[1,1,1,K])).^2./repmat(max(bg,1e-3),[1,1,1,K]),3);
    
    %Warp back into video frames, 
    [h2,w2,~,~] = size(imgs);
    M = zeros([size(imgs,1), size(imgs,2), 1, K]);
    parfor i=1:K
        M(:,:,:,i) = imtransform(moseg(:,:,:,i), maketform('projective', H{i}'), ...
            'XData', [1 w2], 'YData', [1 h2], 'XYScale',1);
    end
    
    %Make binary and fill small holes
    threshold = 0.01*sqrt(h*w);
    parfor i=1:size(moseg,4) 
        M(:,:,:,i) = imfill(medfilt2(M(:,:,:,i),[5,5])>threshold, 'holes');
    end
end


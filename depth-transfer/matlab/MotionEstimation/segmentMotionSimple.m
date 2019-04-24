function moseg = segmentMotionSimple( imgs, flows )
%SEGMENTMOTIONSIMPLE Performs a binary segmentation of moving pixels
% Estimates pixels moving relative to the camera by thresholding the 
% optical flow. 
%
% Input:
%  imgs  - Image/video data of size [height x width x {1,3} x numFrames]
%  flows - Optical flow data corresponding to imgs
% Output:
%  moseg - Binary motion segmentation of estimated moving pixels in imgs of 
%          size [height x width x 1 x numFrames]
%%%%%%%%%%%   Begin segmentMotionSimple   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [h,w,~,K] = size(imgs);
    moseg = sqrt(sum(flows.^2,3));
    
    %Make binary and fill small holes
    threshold = 0.005*sqrt(h*w);
    parfor i=1:K 
        moseg(:,:,:,i) = imfill(medfilt2(moseg(:,:,:,i),[5,5])>threshold, 'holes');
    end
end

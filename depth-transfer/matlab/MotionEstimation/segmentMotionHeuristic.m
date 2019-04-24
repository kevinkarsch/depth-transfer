function moseg = segmentMotionHeuristic( imgs, flows )
%SEGMENTMOTIONHEURISTIC Performs a binary segmentation of moving pixels
% Estimates pixels moving relative to the camera using a heuristic on image
% values and optical flow.
%
% Input:
%  imgs  - Image/video data of size [height x width x {1,3} x numFrames]
%  flows - Optical flow data corresponding to imgs
% Output:
%  moseg - Binary motion segmentation of estimated moving pixels in imgs of 
%          size [height x width x 1 x numFrames]
%%%%%%%%%%%   Begin segmentMotionHeuristic   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imgs = mean(imgs,3);
    [h,w,~,K] = size(imgs);
    
    %Match histograms
    [~,idx] = min(mean(mean(imgs,1),2)); %always match to darkest image
    matchhistto = imhist(imgs(:,:,:,idx));
    B = zeros(h,w,1,K);
    for i=1:K
        B(:,:,:,i) = histeq(imgs(:,:,:,i), matchhistto);
    end
    
    %Create background image by a weighted average of images (based on flow)
    flowweight = 1-sigmoid(sum(flows.^2,3),0.2,100);
    background = sum(imgs.*flowweight,4)./sum(flowweight,4);
    
    %Our motion segmentation heuristic
    moseg = sqrt(sum(flows.^2,3)).*(imgs-repmat(background,[1,1,1,K])).^2 ./ ...
            repmat(background,[1,1,1,K]);

    %Make binary and fill small holes
    threshold = 0.01;
    parfor i=1:size(moseg,4) 
        moseg(:,:,:,i) = imfill(medfilt2(moseg(:,:,:,i),[5,5])>threshold, 'holes');
    end
end


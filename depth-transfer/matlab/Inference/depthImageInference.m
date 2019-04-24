function [depth depth0] = depthImageInference(img, Depth, DepthX, DepthY, ...
                                              SIFTError, Mask, Prior, params)
%DEPTHIMAGEINFERENCE Wrapper for depthVideoInference() for single images
% See depthVideoInference.m for parameter description (parameters are
% identical, except it is assumed that numFrames==1).
%%%%%%%%%%%   Begin depthVideoInference   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [h,w,~,K] = size(img);
    flows = zeros(h,w,2,K);
    FlowError = zeros(size(SIFTError));
    [depth depth0] = depthVideoInference(img, flows, Depth, ...
        DepthX, DepthY, SIFTError, FlowError, Mask, Prior, params);
end


function weights = getTemporalWeights(imgs, temporalGradient, center, scale)
%GETSPATIALWEIGHTS Computes spatial smoothness weights based on an input image or video
% Input:
%  img              - 4D img data [height x width x {1,3} x numFrames]
%  temporalGradient - Temporal gradient operator of size 
%                     [height*width, height*width], usually defined with
%                     [~,~,temporalGradient] = getVideoGradients(flow)
% Output:
%  weights - 4D weights [height x width x 1 x numFrames]; higher implies
%            higher smoothness confidence
%%%%%%%%%%%   Begin getTemporalWeights   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if( ~exist('center','var') )
        center = 0.005;
    end
    if( ~exist('scale','var') )
        scale = 1000;
    end
    
    [h,w,~,K] = size(imgs);
    N = h*w*K;
    B = mean(imgs,3);
    %Linear optical flow operator
    F = -(temporalGradient-spdiags(diag(temporalGradient),0,N,N));
    %Reprojection error from flow
    absFlowErr = abs(imfilter(reshape((F*B(:)-B(:)),[h,w,1,K]), ... 
                     fspecial('gaussian',[5,5],1)));
    weights = 1-sigmoid(absFlowErr, center, scale); %Scale the error nonlinearly
    weights(:,:,:,end) = 0; %Unreliable flow info on boundary (ie last frame)
end

function weights = getSpatialWeights(imgs, center, scale)
%GETSPATIALWEIGHTS Computes spatial smoothness weights based on an input image or video
% Input:
%  img    - 4D video/img data [height x width x {1,3} x numFrames]
%  center - Center of sigmoid (for soft thresholding of gradients)
%  scale  - Scale of sigmoid (for soft thresholding of gradients)
% Output:
%  weights - 4D weights [height x width x 2 x numFrames]; higher implies
%            higher smoothness confidence
%%%%%%%%%%%   Begin getSpatialWeights   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if( ~exist('center','var') )
        center = 0.001;
    end
    if( ~exist('scale','var') )
        scale = 5000;
    end

    B = mean(imgs,3);
    for i=1:size(B,4)
        B(:,:,:,i) = medfilt2(B(:,:,:,i),[5,5]);
    end
    weightsx = 1-sigmoid(abs(imfilter(B,[-1 1 0],'same','replicate')).^2, center, scale);
    weightsx(:,1,:,:) = 0; %Unreliable gradient info on boundary
    weightsy = 1-sigmoid(abs(imfilter(B,[-1 1 0]','same','replicate')).^2,center, scale);
    weightsy(1,:,:,:) = 0; %Unreliable gradient info on boundary
    weights = cat(3,weightsx,weightsy);
end

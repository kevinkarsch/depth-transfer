function [lefts, rights, anaglyphs, disparity] = stereoWarp( imgs, depths, flows, displevels)
%STEREOWARP warps an image sequence (with corresponding flow and depth) to a stereo sequence
  %Input:
    %Assuming [h,w] are the video dimensions and K is the number of frames:
    %imgs is of size [h,w,3,K] - video frames
    %depths of size [h,w,1,K] - depth frames
    %flows of size [h,w,2,K] - optical flow results
    %displevels - maximum amount of disparity between resulting stereo frames (in pixels)
  %Output:
    %lefts/rights - estimated stereo frames (same size as imgs)
    %anaglyphs - stereo frames converted to anaglyph
    %disparity - estimated disparity frames

%% Spatial + temporal smoothness parameters (higher => smoother)
    smoothCoeff_x = 20;
    smoothCoeff_y = 4;
    smoothCoeff_tem = 10;
    
%% Prepare optimization
    fprintf('Preparing optimization...');
    preparetime = tic;
    
    [h,w,~,K] = size(imgs);
    if(~exist('flows','var'))
        flows = zeros(h,w,2,K);
    end
    if( ~exist('displevels','var') )
        displevels = floor(sqrt(h*w)/30);
    end
    N = h*w*K;
    gr = mean(imgs,3);
    grs = imfilter(gr, fspecial('gaussian', [5,5],1));
    eps = 1e-4;
    
    lowprc = prctile(depths(:),1);
    depths(depths<lowprc) = lowprc;
    highprc = prctile(depths(:),99);
    depths(depths>highprc) = highprc;
    
    %One heuristic for converting depth to disparity
    disparity0 = imnormalize(1./(1+imnormalize(depths))).*displevels;
    
    [Gx, Gy, Gt] = getVideoGradients(flows);
    F = -(Gt-spdiags(diag(Gt),0,N,N)); %Linear optical flow operator
    absFlowErr = abs(imfilter(reshape((F*grs(:)-grs(:)),[h,w,1,K]), fspecial('gaussian',[5,5],1)));
    flowConf = 1-sigmoid(absFlowErr,0.005,1000); %Scale the error nonlinearly
    flowConf(:,:,:,end) = 0; %No flow info for last frame
    flowConf(flowConf<eps) = eps;
    wt = spdiags(flowConf(:),0,N,N);
    
    dx = imfilter(grs, [-1 1 0], 'replicate', 'same');
    dy = imfilter(grs, [-1 1 0]', 'replicate', 'same');
    W = (imnormalize(disparity0) + 1.*sigmoid(sqrt(dx.^2+dy.^2),0.01,500))./2;

    A = spdiags(W(:),0,N,N) + smoothCoeff_x.*(Gx'*Gx) + smoothCoeff_y.*(Gy'*Gy) + smoothCoeff_tem.*(Gt'*wt*Gt);
    b = W(:).*disparity0(:);
    
    fprintf('done! (%.5fs)\n', toc(preparetime));
    
%% Optimize
    fprintf('Optimizing warps...\n');
    opttime = tic;
    
    %x = A\b;
    fprintf('\tConstructing preconditioner...\n');
    M = ichol(A,struct('michol','on'));
    fprintf('\tOptimizing...\n\t');
    x = pcg(A,b,1e-6,150,M,M',disparity0(:));
    
    disparity = reshape(x,[h,w,1,K]);
    warpleft = -disparity./2;
    warpright = disparity./2;
    
    fprintf('Done! (%.5fs)\n', toc(opttime));
    
%% Warp and save results
    fprintf('Saving results...');
    savetime = tic;

    lefts = zeros(size(imgs));
    rights = zeros(size(imgs));
    for i=1:K
        lefts(:,:,:,i) = clip(warpImage(imgs(:,:,:,i), warpleft(:,:,:,i), zeros(h,w)));
        rights(:,:,:,i) = clip(warpImage(imgs(:,:,:,i), warpright(:,:,:,i), zeros(h,w)));
    end
    offset = floor(max(disparity(:)).*.8);
    sidebysides = [lefts(:,offset+1:end,:,:), rights(:,1:end-offset,:,:)];
    
    anaglyphs = zeros(size(makeAnaglyph(sidebysides(:,:,:,1))));
    for i=1:K
        anaglyphs(:,:,:,i) = clip(makeAnaglyph(sidebysides(:,:,:,i)));
    end
    
    fprintf('done! (%.5fs)\n\n', toc(savetime));
end

function [ I_norm ] = imnormalize( img, low_prc, up_prc )
	I = img;
    if ~isfloat(I)
        I = double(I);
    end
	if( abs(max(I(:))-min(I(:)))<1e-8 )
        I_norm = img;
    else
        if(nargin==1)
            I_norm = (I - min(I(:))) / (max(I(:)) - min(I(:)));
        else
            m = prctile(I(:),low_prc);
            M = prctile(I(:),up_prc);
            I_norm = (I-m)/(M-m);
        end
	end
end

function sx = sigmoid(x, center, scale)
    sx = 1./(1+exp(scale.*(center-x)));
end

function cv = clip(v, bds)
    if(nargin<2)
        bds = [0,1];
    end
	cv = v;
    cv(v(:)<bds(1)) = bds(1);
    cv(v(:)>bds(2)) = bds(2);
end

function [depth depth0] = depthVideoInference(imgs, flows, motionseg, Depth, DepthX, DepthY, ...
                                              SIFTError, FlowError, Mask, Prior, params)
%DEPTHVIDEOINFERENCE Optimization procedure for Depth Transfer inference
% See depthTransfer.m for a demonstration of setting arguments/parameters
%
% Input:
%  imgs        - 4D video data of size [height x width x {1,3} x numFrames]
%                for which depth should be inferred
%  flows       - 4D optical flow data where flow(:,:,:,i) is the optical 
%                flow from imgs(:,:,:,i) to imgs(:,:,:,i+1). flow(:,:,1,:)
%                encodes the horizontal flow component, and flow(:,:,2,:)
%                encodes the vertical flow component. Of size
%                [height x width x 2 x numFrames], and flow(:,:,:,end)
%                is typically zeros(height,width,2).
%  motionseg   - Binary motion segmentation where motionseg(i)==1 implies
%                the i^{th} pixel is hypothesized to be moving relative to
%                the camera. See 'matlab/MotionSeg' for a few motion
%                estimation algorithms. One way to compute this is with:
%                   > motionseg = segmentMotionHeuristic(imgs,flows);
%  Depth       - Warped candidate depth maps per each frame of input. Of
%                size [[height x width x numCandidates x numFrames].
%                Depth(:,:,i,j) corresponds to the i^{th} warped candidate
%                for frame imgs(:,:,:,j). This warping is done using
%                SIFTflow (see depthTransfer.m), and in pseudocode is:
%                   > imgSIFT = densesift(imgs(:,:,,:,j);
%                   > candidateSIFT = densesift(candidate{i}.img)
%                   > warp = SIFTflow(imgSIFT, candidateSIFT)
%                   > Depth(:,:,i,j) = warp(candidate{i}.depth)
%  DepthX      - Warped candidate depth derivatives (horizontal). Same size
%                as Depth. Pseudocode:
%                   > imgSIFT = densesift(imgs(:,:,,:,j);
%                   > candidateSIFT = densesift(candidate{i}.img)
%                   > warp = SIFTflow(imgSIFT, candidateSIFT)
%                   > depthX = imfilter(candidate{i}.depth, [-1, 1, 0]);
%                   > DepthX(:,:,i,j) = warp(depthX)
%  DepthY      - Warped candidate depth derivatives (vertical). Same size
%                as Depth. Pseudocode:
%                   > imgSIFT = densesift(imgs(:,:,,:,j);
%                   > candidateSIFT = densesift(candidate{i}.img)
%                   > warp = SIFTflow(imgSIFT, candidateSIFT)
%                   > depthY = imfilter(candidate{i}.depth, [-1; 1; 0]);
%                   > DepthY(:,:,i,j) = warp(depthY)
%  SIFTError   - SIFT flow warping error (greater SIFTerror implies greater 
%                warping/SIFTflow errors). Same size as Depth. Can be 
%                computed as reprojection error. Pseudocode:
%                   > imgSIFT = densesift(imgs(:,:,,:,j);
%                   > candidateSIFT = densesift(candidate{i}.img)
%                   > warp = SIFTflow(imgSIFT, candidateSIFT)
%                   > reproj = imgSIFT - warp(candidateSIFT)
%                   > SIFTError(:,:,i,j) = sum( abs(reproj./255), 3)
%  FlowError   - Optical flow warping error. Same size as depth. Can be
%                computed as reprojection error. Pseudocode:
%                   > flow = opticalflow(img(:,:,:,j), img(:,:,:,j+1))
%                   > reproj = img(:,:,:,j) - flow(img(:,:,:,j+1))
%                   > FlowError(:,:,i,j) = sum( abs(reproj), 3)
%  Mask        - Warped candidate depth hole masks. Same size as Depth. 
%                If training depths contain some pixels with unknown or 
%                interpolated depth (e.g. holes in Kinect data), Mask tells 
%                the optimization to not use these pixels (Mask(i)==0 
%                implies Depth{X,Y}(i) are unreliable).  Make sure to warp 
%                the hole masks as with other inputs. Pseudocode:
%                   > imgSIFT = densesift(imgs(:,:,,:,j);
%                   > candidateSIFT = densesift(candidate{i}.img)
%                   > warp = SIFTflow(imgSIFT, candidateSIFT)
%                   > Mask(:,:,i,j) = warp(candidate{i}.holeMask)
%  Prior       - Depth prior used during optimization of size
%                [height x width x 1 x numFrames]. Can be computed as:
%                   > Prior = repmat(computePrior(...), [1,1,1,numFrames]);
%  params(=[]) - Optimization parameters. See initializeProject.m for list
%                and descriptions (or read "%Optimization params" below).
%
% Output:
%  depth  - Inferred depth map corresponding to imgs. Of size
%           [[height x width x 1 x numFrames]
%  depth0 - Initial inferred depth map prior to optimization (for now, this
%           is set as depth0 = median(Depth,3);). Same size as depth.
%
%%%%%%%%%%%   Begin depthVideoInference   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totaltime = tic;
setuptime = tic;
fprintf('\nSetting up inference...');
%Optimization params
try
    alpha = params.alpha;
    gradCoeff = params.gradCoeff;
    smoothCoeff_spa = params.smoothCoeff_spa;
    smoothCoeff_tem = params.smoothCoeff_tem;
    priorCoeff = params.priorCoeff;
    motionCoeff = params.motionCoeff;
catch %#ok<CTCH>
    alpha = 1; %Mixes sift and flow matching error confidence weights (higher => more emphasis on flow error)
    gradCoeff = 10; %Mixes absolute and depth gradients (higher => more ephasis on gradients)
    smoothCoeff_spa = 10.0; %Spatial smoothness coefficient
    smoothCoeff_tem = 100.0; %Temporal smoothness coefficient
    priorCoeff = 0.5; %Larger => enforce depth prior more strongly
    motionCoeff = 5.0; %Larger => enforce motion prior more strongly
end
eps = 1e-3;

%Initialize depth
depth0 = median(Depth,3);
depth = depth0;

%Compute operators/constants
[h,w,~,K] = size(imgs);
N = h*w*K;
M = size(Depth,3); 
[Gx, Gy, Gt] = getVideoGradients( flows ); %x,y, time gradients
matchingWeightSIFT = 1-sigmoid(imnormalize(SIFTError),0.5,10);
matchingWeightFlow = 1-sigmoid(imnormalize(FlowError),0.5,10);
matchingWeight = (matchingWeightSIFT + alpha.*matchingWeightFlow)./2;
spatialWeight = getSpatialWeights(imgs);
flowConf = getTemporalWeights(imgs, Gt);

%Use a direct solver (\) if the system is small(ish)
useDirect = true; 
if(N>500000) 
    useDirect = false;
end

%Permute dimensions for easy reshapes
Depth = permute(Depth, [1,2,4,3]);
DepthX = permute(DepthX, [1,2,4,3]);
DepthY = permute(DepthY, [1,2,4,3]);
Mask = permute(Mask, [1,2,4,3]);
matchingWeight = permute(matchingWeight, [1,2,4,3]);

%Make all weights and costs strictly positive (ensure system is SPD)
depth(depth(:)<eps) = eps;
Depth(Depth(:)<eps) = eps;
Prior(Prior(:)<eps) = eps;
flowConf(flowConf(:)<eps) = eps;
matchingWeight(matchingWeight(:)<eps) = eps;
spatialWeight(spatialWeight(:)<eps) = eps;
Mask(Mask(:)<eps) = eps;
motionseg(motionseg(:)<eps) = eps;

fprintf('done. [%6.2fs]\n', toc(setuptime));

%Optimize
fprintf('Optimizing...\n');
if(useDirect)
    nOuterIterations = 10;
    fprintf('    %14s | %14s\n', 'IRLS iteration','iteration time'); 
else
    nOuterIterations = 10;
    fprintf('    %14s | %12s | %8s | %14s\n', 'IRLS iteration', 'CG residual', 'CG itrs', 'iteration time'); 
end
for i=1:nOuterIterations
    itrtime = tic;
    %Data term
    dw = reshape(matchingWeight,[N,M]);
    [~,W1] = phi( repmat(depth(:),[1,M])-reshape(Depth,[N,M]) );
    W1 = W1 .* dw .* reshape(Mask,[N,M]);
    [~,W1x] = phi( repmat(Gx*depth(:),[1,M])-reshape(DepthX,[N,M]) );
    W1x = W1x .* gradCoeff.*dw .* reshape(Mask,[N,M]);
    [~,W1y] = phi( repmat(Gy*depth(:),[1,M])-reshape(DepthY,[N,M]) );
    W1y = W1y .* gradCoeff.*dw .* reshape(Mask,[N,M]);
    
    %Smoothness
    sx = reshape(spatialWeight(:,:,1,:),[N,1]);
    [~,W2x] = phi( Gx*depth(:) );
    W2x = W2x .* sx.*smoothCoeff_spa;
    sy = reshape(spatialWeight(:,:,2,:),[N,1]);
    [~,W2y] = phi( Gy*depth(:) );
    W2y = W2y .* sy.*smoothCoeff_spa;
    st = flowConf(:);
    [~,W2t] = phi( Gt*depth(:) );
    W2t = W2t .* st.*smoothCoeff_tem;
    
    %Prior
    [~,W3] = phi( depth(:)-Prior(:) );
    W3 = W3 .* priorCoeff;
    
    %Motion prior
    MotionPr = getMotionPrior(depth, motionseg);
    [~,W4] = phi( depth(:) - MotionPr(:) );
    W4 = W4 .* motionseg(:) .* motionCoeff;
    
    %Downweight motion-segmented parts
    if(motionCoeff>0)
        W1 = W1 .* (1-repmat(motionseg(:),[1,M]));
        W1x = W1x .* (1-repmat(motionseg(:),[1,M]));
        W1y = W1y .* (1-repmat(motionseg(:),[1,M]));
        W3 = W3 .* (1-motionseg(:));
    end
    
    %Create system
    A = spdiags(sum(W1,2),0,N,N) + ...
        Gx'*spdiags(sum(W1x,2),0,N,N)*Gx + ...
        Gy'*spdiags(sum(W1y,2),0,N,N)*Gy + ...
        Gx'*spdiags(W2x,0,N,N)*Gx + ...
        Gy'*spdiags(W2y,0,N,N)*Gy + ...
        Gt'*spdiags(W2t,0,N,N)*Gt + ...
        spdiags(W3,0,N,N) + ...
        spdiags(W4,0,N,N);
    b = sum(W1.*reshape(Depth,[N,M]),2) + ...
        Gx'*sum(W1x.*reshape(DepthX,[N,M]),2) + ...
        Gy'*sum(W1y.*reshape(DepthY,[N,M]),2) + ...
        W3.*Prior(:) + ...
        W4.*MotionPr(:);

    %Solve
    if(useDirect)
        depth = reshape(A\b, size(depth));
        fprintf('    %8d    %15.5fs\n', i, toc(itrtime));
    else %Otherwise, use an iterative method
        try
            P = ichol(A, struct('michol','on'));
        catch %#ok<CTCH>
            P = speye(N,N);
            warning(['Cholesky factorization failed, possibly because' ...
                     'the system is not positive definite. This' ...
                     'usually means something very wrong has happened!\n']); %#ok<WNTAG>
        end
        [x,flag,relres,iter,resvec] = pcg(A,b,1e-6,300,P,P',depth(:)); %#ok<ASGLU,NASGU>
        depth = reshape(x, size(depth));
        fprintf('    %8d %19.3e %8d %15.5fs\n', i, relres, iter, toc(itrtime));
    end
end
fprintf('done. [%6.2fs]\n\n', toc(totaltime));
end

function [v, gradv] = phi(x)
    eps = 0.01^2;
    v = sqrt(x.^2+eps);
    gradv = 1./v;
end

function mp = getMotionPrior(depth, motionseg)
    [h,w,~,K] = size(depth);
    Y = ndgrid(1:h,1:w);
    mp = zeros(h,w,1,K);
    for i=1:K
        [L, num] = bwlabel(motionseg(:,:,:,i)>0.5,8);
        for j=1:num
            seg = double(L==j);
            [~,idx] = max(Y(:).*seg(:));
            mp(:,:,:,i) = mp(:,:,:,i) + seg.*repmat(depth(idx+(i-1)*h*w), [h,w,1]);
        end
    end
end

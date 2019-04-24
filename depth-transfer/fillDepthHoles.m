function depth = fillDepthHoles( img, depth0, flow, holeMask )
%FILLDEPTHHOLES Fills holes in depth image(s) using appearance information
%
% Input:
%  img          - 4D video/img data [height x width x {1,3} x numFrames]
%  depth0       - 4D depth data [height x width x 1 x numFrames] (should be
%                 consistent with img)
%  flow         - 4D flow data [height x width x 2 x numFrames]. If
%                 flow = [], temporal info is not used (slower but less
%                 accurate)
%  holeMask     - 4D binary hole data [height x width x 1 x numFrames] 
%                 (should be consistent with img). If a pixel in holeMask
%                 is 0 (or false), it will be interpolated
%
% Output:
%  depth        - filled 4D depth data [height x width x 1 x numFrames]
%                 (holes defined from holeMask will be filled)
%
%%%%%%%%%%%   Begin fillDepthHoles   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Error checking
    [h,w,~,K] = size(img);
    try
        depth0 = reshape(depth0, [h,w,1,K]);
        flow = reshape(flow, [h,w,2,K]);
        holeMask = reshape(holeMask, [h,w,1,K]);
    catch %#ok<CTCH>
        error('All inputs must have the same spatial dimensions and number of frames');
    end
	
    %Initialize operators
    holeMask = (holeMask ~= 0);
    N = w*h*K;
    [Gx, Gy, Gz] = getVideoGradients(flow);
    B = mean(img,3);
    
    spatialWeights = getSpatialWeights(B, 0.01, 500);
    spatialWeights(0<spatialWeights & spatialWeights<0.01) = 0.01;
    Wx = spdiags(reshape(spatialWeights(:,:,1,:),[N,1]), 0,N,N);
    Wy = spdiags(reshape(spatialWeights(:,:,2,:),[N,1]), 0,N,N);
    
    flowConf = getTemporalWeights(B, Gz, 0.001, 5000);
    Wz = spdiags(flowConf(:),0,N,N);

    %Make system
    M = double(spdiags(holeMask(:),0,N,N));
    M(:,sum(M,2)==0) = [];
    d0 = depth0(holeMask(:));
    A = Gx'*Wx*Gx + Gy'*Wy*Gy + Gz'*Wz*Gz;

% Hard constraints => too much memory and PCG requires A to be SPD
%     %Solve constrained quadratic system (lagrangian approach)
%     % argmin_{x} x'*A*x, s.t. Mx = d0
%     % ends up solving: argmin_{x,l} L = x'*A*x + (Mx-d)*l
%     % => L_x = 2*Ax + Ml = 0, L_l = Mx-d = 0
%     % => AL = [2*A M; M 0]*[x;l] = [0;d] = bL
%     AL = [2*A M; M' sparse(size(M,2),size(M,2))];
%     bL = [sparse(N,1); d0];
%     x0 = [depth0(:); sparse(size(M,2),1)];
%     P = speye(size(AL)); %ichol(AL,struct('michol','on'));
%     %x = pcg(AL,bL,1e-6,150,P',P, x0);
%     x = AL\bL;

% Instead use soft constraints
    %Solve argmin_{x} x'*smoothWeight*A*x + sum((Mx-d0)^2)
    smoothWeight = 1e-2;
    P = speye(size(A)); %ichol(smoothWeight.*A + M*M',struct('michol','on'));
    [x,~] = pcg(smoothWeight.*A + M*M', M*d0, 1e-6,300, P', P, depth0(:));

    %Reformat depth to same shape/type as depth0
    depth = cast(reshape(full(x(1:N)),size(depth0)), class(depth0));
end

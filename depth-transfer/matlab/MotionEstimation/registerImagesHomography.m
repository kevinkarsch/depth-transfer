function [I2warp, H] = registerImagesHomography(I1,I2,t)
%REGISTERIMAGESHOMOGRAPHY Registers two images together via homography
% Input:
%  I1       - First image to be registered
%  I2       - Second image to be registered
%  t[=1e-3] - Threshold to discard matches
% Output:
%  I2warp - Image I2 warped to I1's size/coordinate frame (cropped if
%           necessary)
%  H      - The computed homography that warps I2 to I1
%
% NOTE: This function relies on a few of Peter Kovesi's functions for
%  computer vision. For the original functions, please see Peter's site:
%          http://www.csse.uwa.edu.au/~pk/research/matlabfns/
% The necessary functions have been included with this package at the
% default location of 'kovesi_fns' within this directory. If you change the
% path, make sure to modify the variable(s) below.
registerImagesHomography_m_path = fileparts(mfilename('fullpath')); %Path to this function
KOVESI_FNS_PATH = fullfile(registerImagesHomography_m_path, 'kovesi_fns');
%%%%%%%%%%%   Begin registerImagesHomography   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(nargin<3)
        t = 1e-3; % Distance threshold for deciding outliers
    end
    
    %Check if KOVESI_FNS_PATH is already on MATLAB's path
    pathAlreadySet = any(strcmpi(KOVESI_FNS_PATH, regexp(path, pathsep, 'split')));
    if(~pathAlreadySet) %Add it otherwise
        addpath(KOVESI_FNS_PATH);
    end
    
    %Matching harris corners didn't work so well
    %Ig1 = rgb2gray(I1);
    %Ig2 = rgb2gray(I2);
    %[~, r1, c1] = harris(Ig1.*255, 1, 50, 3);
    %[~, r2, c2] = harris(Ig2.*255, 1, 50, 3);
	%[m1,m2] = matchbymonogenicphase(Ig1, [r1';c1'], Ig2, [r2';c2'], 11, 50, 1, 10, 4, 0.2);
    
    %Using SIFT instead (see findMatchingFeatures.m)
    threshold = 1.5;
    [m1,m2] = findMatchingFeatures(I1, I2, threshold); 
    
    x1 = [m1; ones(1,length(m1))];
    x2 = [m2; ones(1,length(m1))];
    if(size(x1,2)>=4)
        H = ransacfithomography(x1, x2, t);
    else
        warning('registerImagesHomography: Not enough matches detected'); %#ok<WNTAG>
        H = eye(3);
    end
	
	I2warp = imtransform(I2, maketform('projective', inv(H')), ...
        'XData', [1 size(I1,2)], 'YData', [1 size(I1,1)], 'XYScale',1);
    
    if(~pathAlreadySet) %Remove from path if necessary
        rmpath(KOVESI_FNS_PATH);
    end
end

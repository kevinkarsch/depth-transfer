function [m1, m2] = findMatchingFeatures(I1, I2, threshold)
%FINDMATCHINGEATURES Computes matching features between two images
% Input:
%  I1              - First image to match
%  I2              - Second image to match
%  threshold(=1.5) - Distance threshold for keeping matching points
% Output:
%  m1 - Matches in the first image
%  m2 - Matches in the second image
%
% NOTE: This code relies on Andrea Vedaldi's SIFT implementation, which is
%  available here: http://www.vlfeat.org/~vedaldi/code/sift.html.
%  Precompiled binaries have been provided with this package, but if you
%  must recompile, see the 'sift' directory located in the same directory
%  as this findMatchingFeatures.m. If you move the sift directory, make
%  sure to modify the below path:
findMatchingFeatures_m_path = fileparts(mfilename('fullpath')); %Path to this function
VEDALDI_SIFT_PATH = fullfile(findMatchingFeatures_m_path, 'sift');
%%%%%%%%%%%   Begin findMatchingFeatures   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(~isfloat(I1))
        I1 = im2double(I1);
    end
    if(~isfloat(I2))
        I2 = im2double(I2);
    end
    if(nargin<3)
        threshold = 1.5;
    end
    
    I1 = rgb2gray(I1);
    I2 = rgb2gray(I2);
    
    %Check if VEDALDI_SIFT_PATH is already on MATLAB's path
    pathAlreadySet = any(strcmpi(VEDALDI_SIFT_PATH, regexp(path, pathsep, 'split')));
    if(~pathAlreadySet) %Add it otherwise
        addpath(VEDALDI_SIFT_PATH);
    end
    
    %Find features using Vedaldi's SIFT
    [framesBase,descriptorsBase] = sift(I1, 'Verbosity', 0);
    [framesInput,descriptorsInput] = sift(I2, 'Verbosity', 0);
    
    %Match feature points
    matches = siftmatch(descriptorsBase, descriptorsInput, threshold);
    m1 = framesBase(1:2, matches(1,:));
    m2 = framesInput(1:2, matches(2,:));

    %Remove from path if necessary
    if(~pathAlreadySet)
        rmpath(VEDALDI_SIFT_PATH);
    end
end


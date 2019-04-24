function gist = im2gist(im)
% This (and the contents of 'gist/') is a modified version of outdated GIST 
%  code corresponding to the paper [Oliva and Torralba IJCV'01]. A newer 
%  version of the code can be found at: 
%       http://people.csail.mit.edu/torralba/code/spatialenvelope/

addpath('gist');

% Default parameters
param.imageSize = 128;
param.orientationsPerScale = [8 8 8 8];
param.numberBlocks = 4;
param.fc_prefilt = 4;

if(nargin==0) %Return number of features if no args given
    gist = param.numberBlocks^2 *sum(param.orientationsPerScale);
    return;
end

param.G = createGabor(param.orientationsPerScale, param.imageSize);

if size(im,3) >1
    im = rgb2gray(im);
end
img = single(im);

% resize and crop image to make it square
img = imresizecrop(img, param.imageSize, 'bilinear');

% scale intensities to be in the range [0 255]
img = img-min(img(:));
img = 255*img/max(img(:));

% prefiltering: local contrast scaling
output    = prefilt(img, param.fc_prefilt);

% get gist:
gist = gistGabor(output, param.numberBlocks, param.G);


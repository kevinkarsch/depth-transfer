% function to demo the usage of mexDenseSIFT.m

% load an RGB image
im = imread('Mars-1.jpg');

% set parameters
% cellSize = [2,4];
cellSize = 3;
stepSize = 1;

% call the mex function to convert to dense SIFT image;
tic; sift = mexDenseSIFT(im,cellSize,stepSize); toc;

% display the SIFT image
figure;imshow(showColorSIFT(sift));

% show the histogram of the magnitude of each SIFT desciptor
% it's not strictly 1 because a normalization is added to avoid 0
foo = double(sift)/255;
foo = sqrt(sum(foo.^2,3));

fprintf('mean: %f, std: %f\n',mean(foo(:)),std(foo(:)));
figure;
hist(foo(:));

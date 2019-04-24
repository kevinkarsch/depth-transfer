%Load data
imgs = cell(1,1,1,3);
depths = cell(1,1,1,3);
for i=1:3
    imgs{i} = im2double(imread(fullfile('data',sprintf('car%d.png',i))));
    tmp = load(fullfile('data',sprintf('car%d-depth.mat',i)));
    depths{i} = tmp.depth;
end

%Compute forward optical flow: flows are 2-channel images, i.e. flow(:,:,1) is horizontal, flow(:,:,2) is vertical
flows = cell(size(imgs));
for i=1:numel(imgs)-1
    flows{i} = opticalflow(imgs{i}, imgs{i+1});
end
flows{end} = zeros(size(flows{1})); %Don't know flow for last frame

%Convert to 4D arrays (as stereoWarp.m expects)
imgs = cell2mat(imgs);
depths = cell2mat(depths);
flows = cell2mat(flows);

%Run stereoWarp
[lefts, rights, anaglyphs] = stereoWarp(imgs, depths, flows, 10);

%Save results
if(~exist('results','dir')); mkdir('results'); end
for i=1:3
    imwrite(lefts(:,:,:,i), fullfile('results',sprintf('left%d.png',i)));
    imwrite(rights(:,:,:,i), fullfile('results',sprintf('right%d.png',i)));
    imwrite(anaglyphs(:,:,:,i), fullfile('results',sprintf('anaglyph%d.png',i)));
end
    

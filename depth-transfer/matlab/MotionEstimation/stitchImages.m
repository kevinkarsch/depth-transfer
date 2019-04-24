function [stitch, img_overlap, mask_overlap] = stitchImages( I )
%STITCHIMAGES Performs panoramic stitching on an image sequence (no bundle adjustment)
% Input:
%  I  - Image/video data of size [height x width x {1,3} x numFrames]
% Output:
%  stitch       - resulting panorama of size [height2 x width2 x {1,3}]
%  img_overlap  - A sequence of all warped/aligned frames. This is the
%                 result of image alignment prior to median filtering (i.e.
%                 stitch = median_filter(img_overlap, 4). Of size
%                 [height2 x width2 x {1,3} x numFrames].
%  mask_overlap - A sequence of masks indicating valid pixels to aid in
%                 median filtering. Takes a value of 1 if img_overlap 
%                 contains useable pixel data. Same size as img_overlap.
%%%%%%%%%%%   Begin stitchImages   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [h,w,d,k] = size(I);
    
    H = cell(k,1);
    H{1} = eye(3,3);
    xdata = zeros(k,2);
    xdata(1,:) = [1,w];
    ydata = zeros(k,2);
    ydata(1,:) = [1,h];
    
    %Find matching homographies/warpings
    for i=2:k
        [~, H{i}] = registerImagesHomography(I(:,:,:,1), I(:,:,:,i));
        [~, xdata(i,:), ydata(i,:)] = imtransform(I(:,:,:,i), maketform('projective', inv(H{i}')), 'XYScale',1);
    end
    
    %Get bounds of warped images
    bounds_x = [min(xdata(:)), max(xdata(:))];
    bounds_y = [min(ydata(:)), max(ydata(:))];
    
    %Stitch
    img_overlap = zeros(ceil(bounds_y(2)-bounds_y(1)+1), ceil(bounds_x(2)-bounds_x(1)+1),d,k);
    mask_overlap = zeros(size(img_overlap));
    for i=1:k
        img_overlap(:,:,:,i) = imtransform(I(:,:,:,i), maketform('projective', inv(H{i}')), ...
            'XData', bounds_x, 'YData', bounds_y, 'XYScale',1);
        mask_overlap(:,:,:,i) = imtransform(ones(h,w,d), maketform('projective', inv(H{i}')), ...
            'XData', bounds_x, 'YData', bounds_y, 'XYScale',1);
    end
    
    %Average warped frames
    %stitch = sum(img_overlap,4)./max(sum(mask_overlap,4),1); 
    %Median filter instead?
    [h2,w2,d2,~] = size(img_overlap);
    sort_img_overlap = sort(img_overlap,4, 'descend');
    medianidx = max(ceil(sum(mask_overlap,4)/2),1);
    [X,Y,Z] = ndgrid(0:h2-1, 0:w2-1, 0:d2-1);
    stitch = sort_img_overlap(X+Y*h2+Z*h2*w2 + h2*w2*d2.*(medianidx-1) + 1);
end


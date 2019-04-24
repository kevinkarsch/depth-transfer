function fgist = flow2flowgist( flow )
%FLOW2FLOWGIST Computes "gist" of an optical flow frame
% flow is a height x width x 2 matrix of flow vectors; flow(:,:,1) is the 
% horizontal flow, flow(:,:,2) is the vertical flow
    nblocks = [4,4]; %number of image blocks (uniform grid)
    ftsperblock = 8; %number of features per block
    fgist = zeros(nblocks(1)*nblocks(2)*ftsperblock,1);
    if(nargin==0) %if no args, return total # features
        fgist = numel(fgist);
        return;
    end
    [h,w,~] = size(flow);
    step = floor([h,w]./nblocks);
    for i=1:nblocks(1)
        for j=1:nblocks(2)
            flowblock = flow((i-1)*step(1)+(1:step(1)),(j-1)*step(2)+(1:step(2)),:);
            flowblock = reshape(flowblock,[step(1)*step(2),2]);
            fts = [mean(flowblock), std(flowblock), mean(flowblock.^2), std(flowblock.^2)];
            fgist(((i-1)*nblocks(2)+j-1)*ftsperblock+(1:ftsperblock),:) = fts;
        end
    end
end

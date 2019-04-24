function depthPrior = computePrior( project, trainClips )
%COMPUTEPRIOR Computes a depth prior for DepthTransfer given training files
%
% Input:
%  project    - Project info struct created with initializeProject(...)
%  trainClips - Struct array such that trainClips(i).name contains the 
%               i^{th} training data directory (relative to
%               project.path.data). MATLAB's built in dir() function can 
%               create these structs; ex:
%                   trainClips = ...
%                       dir(fullfile(project.path.data, [BASE_NAME '*']);
%               where BASE_NAME is a name shared by all training 
%               directories
%
% Output:
%  depthPrior - Depth prior result of size [project.h, project.w]
%
%%%%%%%%%%%   Begin computePrior   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Load all depth files from trainClips, average across each clip
    depth = zeros(project.h, project.w, numel(trainClips));
    mask = zeros(project.h, project.w, numel(trainClips));
    for i=1:numel(trainClips)
        trainName = fullfile(project.path.data, trainClips(i).name);
        [~, depthi, ~, maski] = loadData(trainName, [], [project.h, project.w]);
        if(iscell(depthi))
            depthi = cell2mat(reshape(depthi,1,1,1,[]));
            maski = cell2mat(reshape(maski,1,1,1,[]));
        end
        %Weighted average per video (only use valid depth pixels)
        mask(:,:,i) = any(maski,4);
        depth(:,:,i) = sum(double(maski).*depthi,4)./max(sum(double(maski),4),1);
    end
    %Then average across the entire training set
    depthPrior = sum(double(mask).*depth,3)./max(sum(double(mask),3),1);
    %Fill in any holes (very unlikely that the prior contains any though)
    depthPrior(~any(mask,3)) = mean(depth(:));
end


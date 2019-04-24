%MAKE Compiles (mex) necessary files for SIFTflow
function make()
    %This package comes precompiled for many systems, and might not require
    % any mexing. However, if you do need to recompile, make sure your mex
    % environment is setup by running: 'mex -setup'
    %If compile fails:
    %  1) Make sure SIFT_FLOW_DIR is set properly
    %  2) Read the READMEs in each folder of matlab/SIFTflow

    curdir = pwd();
    SIFT_FLOW_DIR = 'matlab/SIFTflow';
    
    try
        fprintf('Compiling mexWarpImageInt...');
        cd( fullfile(SIFT_FLOW_DIR, 'warpImageInt') );
        mex mexWarpImageInt.cpp
        fprintf('done.\n');
    catch %#ok<CTCH>
        fprintf('FAILED (check make.m).\n');
    end
    cd(curdir);
    
    try
        fprintf('Compiling mexDenseSIFT...');
        cd( fullfile(SIFT_FLOW_DIR, 'mexDenseSIFT') );
        mex mexDenseSIFT.cpp Matrix.cpp Vector.cpp
        fprintf('done.\n');
    catch %#ok<CTCH>
        fprintf('FAILED (check make.m).\n');
    end
    cd(curdir);
    
    try
        fprintf('Compiling mexDiscreteFlow...');
        cd( fullfile(SIFT_FLOW_DIR, 'mexDiscreteFlow') );
        mex mexDiscreteFlow.cpp BPFlow.cpp Stochastic.cpp
        fprintf('done.\n');
    catch %#ok<CTCH>
        fprintf('FAILED (check make.m).\n');
    end
    cd(curdir);
    
end

function add_training_data(project)
%ADD_TRAINING_DATA Helper function for adding the demo/example training data

%Reads in a few of the example img/depth pairs and adds them to the
% project's data directory 
trainFilesDir = fullfile('sample_training_data');
trainFiles = dir(fullfile(trainFilesDir, 'img-*.jpg'));
parfor i=1:numel(trainFiles)
    tmpProject = project; %Avoid parfor slicing issue
    [~, name, ~] = fileparts(fullfile(trainFilesDir, trainFiles(i).name));
    basename = name(5:end); %Remove 'img-' prefix
    dataDirName = fullfile(tmpProject.path.data,['Make3D-Train-' basename]);
    if( exist(fullfile(dataDirName, '001'), 'dir') )
        continue; %Training data already exists
    end
    img = imread(fullfile(trainFilesDir, trainFiles(i).name));
    depthFile = dir( fullfile(trainFilesDir, ['depth_sph_corr-' basename '.mat']) );
    foo = load( fullfile(trainFilesDir, depthFile(1).name) );
    depth = foo.Position3DGrid(:,:,4); %Load only depth from laser data
    createData(dataDirName, img, depth, [], false); %false => verbose off
end

%% SET ALL DEFAULT OPTIONS HERE
% check out the README file for detailed instructions (and extra options)

% cd('C:\Users\Shin\Documents\GitHub\Suite2P\local') % start this code in the directory with make_db

% make_db_example;
% make_db
make_db4Alice;
% make_dbM1M2;

[toolbox_path,~,~] = fileparts(which('run_pipeline'));

ops0.toolbox_path = toolbox_path;
if exist(ops0.toolbox_path, 'dir')
	addpath(ops0.toolbox_path) % add local path to the toolbox
else
	error('toolbox_path does not exist, please change toolbox_path');
end

ops0.clustModel  = 'neuropil'; % standard or neuropil
ops0.neuropilSub = 'surround'; % none, surround or model

ops0.useGPU                 = 1; % if you can use a GPU in matlab this accelerate registration approx 3 times
ops0.doRegistration         = 0;
% root paths for files and temporary storage (ideally an SSD drive. my SSD is C)
ops0.RegFileTiffLocation    = []; %'D:/DATA/'; % leave empty to NOT save registered tiffs
ops0.RegFileRoot            = 'E:\Imaging\Suite2P';

ops0.getROIs                = 1;
ops0.getSVDcomps            = 1;
ops0.nSVD                   = 1000; % how many SVD components to keep

ops0.temp_tiff              = '\\research.files.med.harvard.edu\neurobio\HarveyLab\Shin\ShinDataAll\Suite2P\temp\temp.tiff';
ops0.ResultsSavePath        = '\\research.files.med.harvard.edu\neurobio\HarveyLab\Shin\ShinDataAll\Suite2P';

% ops0.temp_tiff              = 'C:\Users\Shin\Documents\MATLAB\ShinDataAll\Suite2P\temp\temp.tiff'; % copy data locally first
% ops0.ResultsSavePath        = 'C:\Users\Shin\Documents\MATLAB\ShinDataAll\Suite2P';
ops0.PhaseCorrelation       = 1; % set to 0 for non-whitened cross-correlation
ops0.SubPixel               = Inf; % 2 is alignment by 0.5 pixel, Inf is the exact number from phase correlation

ops0.showTargetRegistration = 1;
ops0.RootStorage            = '\\research.files.med.harvard.edu\Neurobio\HarveyLab\Shin\ShinDataAll\Imaging'; % location of saved image data
ops0.ShowCellMap            = 1;
ops0.DeleteBin              = 1; % set to 1 to delete bin file after processing

ops0.diameter               = 10; % mean diameter of your ROIs (varies with FOV magnification)

% the following settings shouldn't need to be adjusted
ops0.NavgFramesSVD          = 5000; % how many (pooled) frames to do the SVD based on
ops0.Nk0                    = 1300;  % how many clusters to start with
ops0.Nk                     = 650;  % how many clusters to end with
ops0.nSVDforROI             = 1000;
ops0.niterclustering        = 50;    % how many iterations of clustering

ops0.NimgFirstRegistration  = 500;
ops0.RegPrecision           = 'int16';
ops0.RawPrecision           = 'int16';
ops0.NiterPrealign          = 10;

ops0.LoadRegMean   			= 0; %
ops0.nimgbegend             = 250; % how many frames to average at the beginning and end of each experiment
ops0.sig                    = 0;  % spatial smoothing constant
% (encourages colocalized clusters) OBSOLETE

ops1 = ops0;

% parameters for signal and neuropil calculation

% NOTE: to save time, first run detected ROIs through the gui to select
% "good ROIs", then extract signals and neuropil only on those ROIs;
% to do so, set .getSignal=0 and .getNeuropil=0 in the master_file, then
% run get_signals_and_neuropi with the following parameters
% ops1.processed  = 0; % to load processed file
% ops1.newFile    = 1; % to create new file ending _new.mat, otherwise
%                        processed file is overwritten
% ops1.useSVD     = 1; % has to be SVD as bin-file is only temporary
% set other parameters accordingly


%%
for iexp = 1:length(db)
    for iplane = 1:db(iexp).nplanes
        data_file = fullfile(ops0.ResultsSavePath,db(iexp).mouse_name,db(iexp).date,sprintf('F_%s_%s_plane%d_Nk%d.mat', db(iexp).mouse_name, db(iexp).date, iplane, ops0.Nk));
        if ~exist(data_file,'file')
            
            disp(db(iexp).comments);
            
            % copy files from zserver
            run_pipeline(db(iexp), ops0);
            
            % deconvolved data into (dat.)cl.dcell, and neuropil subtraction coef
            % add_deconvolution(ops1, db(iexp), clustrules);
        end
    end
end
%%

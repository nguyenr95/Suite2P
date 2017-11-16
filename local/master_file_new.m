%% SET ALL DEFAULT OPTIONS HERE

% UPDATE Christmas 2016: number of clusters determined automatically, but
% do specify the "diameter" of an average cell for best results. You can do this with either
% db(iexp).diameter, or ops0.diameter

% check out the README file for detailed instructions (and extra options)
% addpath('D:\CODE\MariusBox\runSuite2P') % add the path to your make_db file

% overwrite any of these default options in your make_db file for individual experiments
make_db; % RUN YOUR OWN MAKE_DB SCRIPT TO RUN HERE

% ops0.toolbox_path = 'C:\CODE\GitHub\Suite2P';
[ops0.toolbox_path,~,~] = fileparts(which('run_pipeline'));

if exist(ops0.toolbox_path, 'dir')
	addpath(genpath(ops0.toolbox_path)) % add local path to the toolbox
else
	error('toolbox_path does not exist, please change toolbox_path');
end

% mex -largeArrayDims SpikeDetection/deconvL0.c (or .cpp) % MAKE SURE YOU COMPILE THIS FIRST FOR DECONVOLUTION

ops0.useGPU                 = 1; % if you can use an Nvidia GPU in matlab this accelerates registration approx 3 times. You only need the Nvidia drivers installed (not CUDA).
ops0.fig                    = 1; % turn off figure generation with 0

% root paths for files and temporary storage (ideally an SSD drive. my SSD is C:/)
ops0.RootStorage            = '\\research.files.med.harvard.edu\Neurobio\HarveyLab\Tier2\Shin\ShinDataAll\Imaging'; % location of saved image data (folder structure, check out README file)
ops0.temp_tiff              = '\\research.files.med.harvard.edu\neurobio\HarveyLab\Tier1\Shin\ShinDataAll\Suite2P\temp\temp.tiff'; % copies each remote tiff locally first, into this file
ops0.RegFileRoot            = 'E:\Imaging\Suite2P';  % location for binary file
ops0.DeleteBin              = 0; % set to 1 for batch processing on a limited hard drive
ops0.ResultsSavePath        = '\\research.files.med.harvard.edu\neurobio\HarveyLab\Tier1\Shin\ShinDataAll\Suite2P'; % a folder structure is created inside
ops0.RegFileTiffLocation    = []; %'D:/DATA/'; % leave empty to NOT save registered tiffs (slow)

% registration options
ops0.doRegistration         = 0; % skip (0) if data is already registered
ops0.showTargetRegistration = 1; % shows the image targets for all planes to be registered
ops0.PhaseCorrelation       = 1; % set to 0 for non-whitened cross-correlation
ops0.SubPixel               = Inf; % 2 is alignment by 0.5 pixel, Inf is the exact number from phase correlation
ops0.NimgFirstRegistration  = 500; % number of images to include in the first registration pass 
ops0.nimgbegend             = 250; % frames to average at beginning and end of blocks

% cell detection options
ops0.ShowCellMap            = 1; % during optimization, show a figure of the clusters
ops0.sig                    = 0.5;  % spatial smoothing length in pixels; encourages localized clusters
ops0.nSVDforROI             = 1000; % how many SVD components for cell clustering
ops0.NavgFramesSVD          = 5000; % how many (binned) timepoints to do the SVD based on
ops0.diameter               = 10; % (added by SK) expected diameter of cells (used for 0.25 * pi/4*diam^2 < npixels < 10*pi/4*diam^2)

% spike deconvolution options
ops0.imageRate              = 30;   % imaging rate (cumulative over planes!). Approximate, for initialization of deconvolution kernel.
ops0.sensorTau              = 2; % decay half-life (or timescale). Approximate, for initialization of deconvolution kernel.
ops0.maxNeurop              = Inf; % for the neuropil contamination to be less than this (sometimes good, i.e. for interneurons)
ops0.recomputeKernel        = 1; % whether to re-estimate kernel during optimization (default kernel is "reasonable", if you give good timescales)
ops0.sameKernel             = 1; % whether the same kernel should be estimated for all neurons (robust, only set to 0 if SNR is high and recordings are long)

% red channel options
% redratio = red pixels inside / red pixels outside
% redcell = redratio > mean(redratio) + redthres*std(redratio)
% notred = redratio < mean(redratio) + redmax*std(redratio)
ops0.redthres               = 1.5; % the higher the thres the less red cells
ops0.redmax                 = 1; % the higher the max the more NON-red cells

%% RUN THE PIPELINE HERE
db0 = db;

for iexp = 1:length(db) %[3:length(db) 1:2]
    for iplane = 1:db(iexp).nplanes
        data_file = fullfile(ops0.ResultsSavePath,db(iexp).mouse_name,db(iexp).date,sprintf('F_%s_%s_plane%d.mat', db(iexp).mouse_name, db(iexp).date, iplane));
        if 1%~exist(data_file,'file')
            
            disp(db(iexp).comments);
            
            run_pipeline(db(iexp), ops0);

            % deconvolved data into (dat.)cl.dcell, and neuropil subtraction coef
            % add_deconvolution(ops0, db0(iexp));
        end
    end
end

%% STRUCTURE OF RESULTS FILE
% cell traces are in dat.Fcell
% neuropil traces are in dat.FcellNeu
% manual, GUI overwritten "iscell" labels are in dat.cl.iscell

% stat(icell) contains all other information
% autoamted iscell label, based on anatomy
% neuropil subtraction coefficient 
% st are the deconvolved spike times (in frames)
% c  are the deconvolved amplitudes
% kernel is the estimated kernel
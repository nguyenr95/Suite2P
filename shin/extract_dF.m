function [dF, traces, rawF, roi, traceNeuropil] = extract_dF(opt,iplane)

if ~isfield(opt,'zoomMicro')
    opt.zoomMicro=2; %fixed zoomMicro
end
if ~isfield(opt,'inNeurop')
    opt.inNeurop=3; %fixed inner diameter of the neuropil mask donut
end
if ~isfield(opt,'outNeurop')
    opt.outNeurop=45; %radius of Neuropil fixed at 45um
end
if ~isfield(opt,'microID')
    opt.microID='b'; %microscope identity
end
if 1 % ~isfield(opt, 'processed') 
    opt.processed = 1; % use processed data in F_...._proc (generated in gui2P)
end
if ~isfield(opt, 'newFile') % save new file '<name>_new.mat', otherwise existing file is overwritten
    opt.newFile = 0;
end
if nargin<2
    if ~isfield(opt,'iplane')
        error('provide iplane')
    else
        iplane = opt.iplane;
    end
end

filenames = dir(sprintf('%sF_%s_%s_plane%d_Nk*.mat',...
    opt.ResultsSavePath, opt.mouse_name, opt.date, iplane));
filenames = {filenames.name};
if isfield(opt, 'Nk')
    ind = cellfun(@strfind, filenames, ...
        repmat({['_Nk' num2str(opt.Nk)]}, size(filenames)), ...
        'UniformOutput', false);
    ind = ~cellfun(@isempty, ind);
    filenames = filenames(ind);
end
ind = cellfun(@strfind, filenames, repmat({'_proc'}, size(filenames)), ...
    'UniformOutput', false);
ind = ~cellfun(@isempty, ind);
if opt.processed == 1
    filenames = filenames(ind);
else
    filenames = filenames(~ind);
end
if length(filenames) > 1
    fprintf(['WARNING: several files for dataset exist\n' ...
        '%s <- using\n'], filenames{1})
    fprintf('%s\n', filenames{2:end})
elseif isempty(filenames)
    error('Could not find cell detection file \n')
end

data = load(fullfile(opt.ResultsSavePath, filenames{1}));
if opt.processed == 1
    data = data.dat;
end
%%
nRoi = sum([data.cl.iscell]);
iRoi = find([data.cl.iscell]);
for r = 1:nRoi
    fprintf('Extracting ROI %03.0f of %03.0f\n', r, nRoi);
    rawF(r,:) = data.F{iplane}.Fcell(iRoi(r),:);
    
    if isfield(roi(r),'indNeuropil') && ~isempty(roi(r).indNeuropil)
        subCoef = roi(r).subCoef;
        indNeuropil = obj.mat2binInd(roi(r).indNeuropil);
        traceNeuropil(r,:) = mean(mov(:, indNeuropil), 2)';
        traces(r,:) = traceCell - traceNeuropil(r,:)*subCoef;
    else
        traces(r,:) = rawF(r,:);
    end
end

dF = dFcalc(traces,rawF,'custom_wfun');
clear mov
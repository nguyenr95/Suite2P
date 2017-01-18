function add_deconvolution(ops, db)
ops = build_ops3(db, ops);
ops0 = ops;

try
    ppool = gcp ;
catch
end

% warning('ops0.imageRate now represents the TOTAL frame rate of the recording over all planes. This warning will be disabled in a future version. ')

for i = 1:length(ops.planesToProcess)
    iplane  = ops.planesToProcess(i);
    
<<<<<<< HEAD
    fpath = sprintf('%s\\F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, iplane, ops.Nk);
    dat = load(fpath);
    
    % fpath = sprintf('%s\\F_%s_%s_plane%d_Nk%d_proc.mat', ops.ResultsSavePath, ...
    %     ops.mouse_name, ops.date, iplane, ops.Nk);
    
%     if exist(fpath, 'file')
%         load(fpath);
%     else
%         fpath = sprintf('%s\\F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ...
%             ops.mouse_name, ops.date, iplane, ops.Nk);
%         dat = load(fpath);
%     end
=======
    fpath = sprintf('%s/F_%s_%s_plane%d_proc.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, iplane);
    if exist(fpath, 'file')
        load(fpath);
    else
        fpath = sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
            ops.mouse_name, ops.date, iplane);
        dat = load(fpath);
    end
>>>>>>> refs/remotes/cortex-lab/master
    
    if isfield(dat, 'dat')
        dat = dat.dat; % just in case...
    end
    
    % overwrite fields of ops with those saved to file
    ops = addfields(ops, dat.ops);
    
    % set up options for deconvolution
    ops.imageRate    = getOr(ops0, {'imageRate'}, 30); % total image rate (over all planes)
    if str2double(db.depth) < 300
        % superficial layers
        ops.sensorTau    = getOr(ops0, {'sensorTau'}, 2); % approximate timescale in seconds
    else
        % deep layers
        ops.sensorTau    = getOr(ops0, {'sensorTau'}, 8); % approximate timescale in seconds
    end
    ops.sameKernel   = getOr(ops0, {'sameKernel'}, 1); % 1 for same kernel per plane, 0 for individual kernels (not recommended)
    ops.sameKernel   = getOr(ops0, {'sameKernel'}, 1);
    ops.maxNeurop    = getOr(ops0, {'maxNeurop'}, Inf);
    ops.recomputeKernel    = getOr(ops0, {'recomputeKernel'}, 1);
    
    
    fprintf('Spike deconvolution, plane %d... \n', iplane)
    
    % split data into batches
    stat = run_deconvolution3(ops, dat);
    
    dat.stat = stat;
    
    fpath = sprintf('%s\\F_%s_%s_plane%d_Nk%d_proc.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, iplane, ops.Nk);
    
    save(fpath, '-struct', 'dat')
end
%

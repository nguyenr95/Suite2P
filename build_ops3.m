function ops = build_ops3(db, ops)

ops.nplanes = getOr(ops, 'nplanes', 1);
ops.nchannels = getOr(ops, 'nchannels', 1);

% ops = db;
if ~iscell(db.mouse_name) 
    % this is the usual case where we have a simple single session recording
    ops = addfields(ops, db);
    
    for k = 1:length(db.expts)
        ops.SubDirs{k}    = num2str(db.expts(k));
    end
    if isempty(db.expts)
        ops.SubDirs{1} = [];
    end
    
    if ~isfield(ops, 'RootDir')
        ops.RootDir = fullfile(ops.RootStorage, ops.mouse_name, ops.date);
        if ops.doRegistration
            ops.RootDir = fullfile(ops.RootStorage, ops.mouse_name, ops.date);
        else
            ops.RootDir = fullfile(ops.RootStorage, ops.mouse_name, ops.date,'Corrected');
            ops.RootDirRaw = fullfile(ops.RootStorage, ops.mouse_name, ops.date);
        end
        disp(ops.RootDir);
    end
    
    % build file list
    ops.fsroot = [];
    for j = 1:ops.nplanes %length(ops.SubDirs) % changed on 16/11/25 by SK
        ops.fsroot{j} = dir(fullfile(ops.RootDir,sprintf('*Slice%02d*.tif',j)));
        if j==1
            ops.rawMovies = dir(fullfile(ops.RootDirRaw,'*tif'));
            for k = 1:length(ops.rawMovies)
                % Do not include overview image
                if isfinite(strfind(ops.rawMovies(k).name,'overview'));
                    overview_ind = k;
                end
            end
            pick_ind = (1:length(ops.rawMovies))~=overview_ind;
            ops.rawMovies = ops.rawMovies(pick_ind);
        end
        % ops.fsroot{j} = dir(fullfile(ops.RootDir,sprintf('*Slice%02d*_File001.tif',j)));
        for k = 1:length(ops.fsroot{j})
            % Do not include overview image
            if isempty(strfind(ops.fsroot{j}(k).name,'overview'));
                ops.fsroot{j}(k).name = fullfile(ops.RootDir, ops.fsroot{j}(k).name);
            else
                ops.fsroot{j}(k).name = [];
            end
        end
        % ops.fsroot{j} = ops.fsroot{j}(1:20); % select some fraction of files
    end
    
    if isfield(db, 'expred') && ~isempty(db.expred) && ...
            (~isfield(db, 'nchannels_red') || isempty(db.nchannels_red))            
        ops.fsred = dir(fullfile(ops.RootDir, num2str(db.expred), '*.tif'));
        for k = 1:length(ops.fsroot{j})
            ops.fsred(k).name = fullfile(ops.RootDir, num2str(db.expred), ops.fsred(k).name);
        end
    end
else
    % here we might have multiple sessions, which we want to be analyzed
    % together (exactly the same FOV)
    nSessions = length(db.mouse_name);
    % a backwards compatible version of db
    dbCompat = db;
    dbCompat.mouse_name = db.mouse_name{1};
    dbCompat.date = db.date{1};
    dbCompat.expts = cell2mat(db.expts(:)');
    ops = addfields(ops, dbCompat);
    ops.db_orig = db;
    
    ops.fsroot = cell(0);
    ops.SubDirs = cell(0);
    for iSession = 1:nSessions
        ops.RootDir = fullfile(ops.RootStorage, db.mouse_name{iSession}, db.date{iSession});
        for iExp = 1:length(db.expts{iSession})
            ops.SubDirs{end+1} = num2str(db.expts{iSession}(iExp));
            ops.fsroot{end+1} = dir(fullfile(ops.RootDir, ops.SubDirs{end}, '*.tif'));
            for iFile = 1:length(ops.fsroot{end})
                ops.fsroot{end}(iFile).name = fullfile(ops.RootDir, ops.SubDirs{end}, ops.fsroot{end}(iFile).name);
            end
        end
    end
    % this line to be backward compatible (just in case)
    ops.RootDir = fullfile(ops.RootStorage, ops.mouse_name, ops.date);
end

try
    % MK code for automatically determining number of planes and channels
    [~, header] = loadFramesBuff(ops.fsroot{1}(1).name, 1, 1, 1);
    
    hh=header{1};
    str = hh(strfind(hh, 'channelsSave = '):end);
    ind = strfind(str, 'scanimage');
    ch = str2num(str(16 : ind(1)-1));
    ops.nchannels = length(ch);
    
    % fastZEnable = sscanf(hh(findstr(hh, 'fastZEnable = '):end), 'fastZEnable = %d');
    fastZEnable = eval(sscanf(hh(findstr(hh, 'fastZEnable = '):end), 'fastZEnable = %s'));
    fastZDiscardFlybackFrames = sscanf(hh(findstr(hh, 'fastZDiscardFlybackFrames = '):end), 'fastZDiscardFlybackFrames = %d');
    if isempty(fastZDiscardFlybackFrames)
        fastZDiscardFlybackFrames = 0;
    end
    stackNumSlices = sscanf(hh(findstr(hh, 'stackNumSlices = '):end), 'stackNumSlices = %d');
    
    ops.nplanes = 1;
    if fastZEnable
        ops.nplanes = stackNumSlices+fastZDiscardFlybackFrames;
    end
    
    str = hh(strfind(hh, 'scanZoomFactor = '):end);
    ind = strfind(str, 'scanimage');
    ops.zoomMicro = str2double(str(18 : ind(1)-1));
    
    % get number of channels of red experiment
    if isfield(db, 'expred') && ~isempty(db.expred) && ...
            (~isfield(db, 'nchannels_red') || isempty(db.nchannels_red))
        [~, header] = loadFramesBuff(ops.fsred(1).name, 1, 1, 1);
        hh=header{1};
        str = hh(strfind(hh, 'channelsSave = '):end);
        ind = strfind(str, 'scanimage');
        ch = str2num(str(16 : ind(1)-1));
        ops.nchannels_red = length(ch);
    end
catch
end

if ~(isfield(ops, 'planesToProcess') && ~isempty(ops.planesToProcess))
    % ops.planesToProcess = 1:ops.nplanes;
    ops.planesToProcess = 1:stackNumSlices; % changed by SK 16/12/14
else
    % planesToProcess is not working right now
    ops.planesToProcess = 1:ops.nplanes;
end

CharSubDirs = '';
for i = 1:length(ops.SubDirs)
    CharSubDirs = [CharSubDirs ops.SubDirs{i} '_'];
end
CharSubDirs = CharSubDirs(1:end-1);
ops.CharSubDirs = CharSubDirs;

ops.ResultsSavePath = sprintf('%s\\%s\\%s\\', ops.ResultsSavePath, ops.mouse_name, ops.date);
% ops.ResultsSavePath = sprintf('%s\\%s\\%s\\%s\\', ops.ResultsSavePath, ops.mouse_name, ops.date, CharSubDirs);


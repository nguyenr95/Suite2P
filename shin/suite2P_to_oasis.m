function suite2P_to_oasis(mouseNum,date_num,varargin)
    % Extract Ca signal from suite2P and deconvolve with constrained_oasisAR1
    
    if ~exist('nslices','var') && ~exist('nplanes','var')
        nSlices = 1;
        sliceSet = 1;
    end
    
    if ~exist('smoothWindow','var') || isempty(smoothWindow)
        smoothWindow = 15;
    end
    
    initials = getInitials(mouseNum);
    
    temp = userpath;
    matlab_dir = temp(1:end-1);
    mouseID = sprintf('%s%03d',initials,mouseNum);
    
    if ~exist('data_file','var')
        folder_name = sprintf(['\\\\research.files.med.harvard.edu\\Neurobio\\HarveyLab\\',...
            'Shin\\ShinDataAll\\Suite2P\\%s\\%d'],mouseID,date_num);
        try
            ops_file = sprintf('regops_%s_%d',mouseID,date_num);
            temp = load(fullfile(folder_name,ops_file));
        catch
            [ops_file,folder_name] = uigetfile([folder_name,'regops*.mat'],'MultiSelect','off');
            temp = load(fullfile(folder_name,ops_file));
        end
        ops0 = temp.ops1{1};
        data_file = dir(fullfile(folder_name,'*_proc.mat'));
    end
    
    movFile = fullfile(ops0.RootDirRaw,'FOV1_001_001.tif'); % SI4?
    if ~exist(movFile,'file')
        movFile = fullfile(ops0.RootDirRaw,'FOV1_00001_00001.tif'); % SI2016
    end
    
    framePeriod = 0.0361;
    
    if ~exist('framePeriod','var')
        castType = 'uint16';
        if 1
            [~, metaDataSI] = tiffRead(movFile,castType);
            if isfield(metaDataSI,'SI4')
                framePeriod = metaDataSI.SI4.scanFramePeriod;
            elseif isfield(metaDataSI,'SI5')
                framePeriod = metaDataSI.SI5.scanFramePeriod;
            elseif isfield(metaDataSI,'SI')
                framePeriod = metaDataSI.SI.hRoiManager.scanFramePeriod;
            else
                warning('Unable to Automatically determine scanFramePeriod')
                framePeriod = input('Input scanFramePeriod: ');
            end
        else
            % for image from Odin rig
            info = imfinfo(movFile);
            temp = info.Software;
            ind = strfind(temp,'SI.');
            temp(ind(2:end)-1) = ';';
            temp(end) = ';';
            eval(temp);
            framePeriod = SI.hRoiManager.scanFramePeriod;
        end
    end
    
    fR = 1/framePeriod;
    tstart = tic;

    for si = 1:nSlices
        if ismember(si,sliceSet)
            
            if nSlices==1
                load(fullfile(folder_name,data_file.name));
            else
                load(fullfile(folder_name,data_file(si).name));
            end
            telapsed = toc(tstart);
            fprintf('loaded slice %02d in %d sec\n',si,round(telapsed))

            if exist('dat','var') % for older files
                img0 = dat.mimg(:,:,5);
                res = dat.res;
                Fcell = dat.Fcell{1};
                FcellNeu = dat.FcellNeu{1};
                isgood = [dat.stat.iscell];
            else
                Fcell = Fcell{1};
                FcellNeu = FcellNeu{1};
            end
            n_cell = length(length(dat.stat));
            pick_cell = find([dat.stat.iscell]);
                        
            startCellInd = 1;
            
            for ci = startCellInd:size(Fcell,1)
                if isgood(ci)
                    %%
                    fprintf('Analyzing ROI %d (%d / %d cells)\n',ci,sum(isgood(1:ci)),sum(isgood));
                    fprintf('Computing Neuropil coefficient ... ')
                    if 0
                        % compute dF
                        dF = Fcell(ci,:) - bsxfun(@times, FcellNeu(ci,:), dat.stat(ci).neuropilCoefficient); % using the method in extractSignals.m
                    else
                        % compute dF/F
                        fit_mode = 'prctile'; %'exp_linear';
                        rawTrace = Fcell(ci,:);
                        subTrace = Fcell(ci,:) - bsxfun(@times, FcellNeu(ci,:), dat.stat(ci).neuropilCoefficient); % using the method in extractSignals.m
                        rawTraceBase = getF_(rawTrace,fit_mode,ops0.imageRate);
                        subTraceBase = getF_(subTrace,fit_mode,ops0.imageRate);
                        dF = (subTrace - subTraceBase)./rawTraceBase;
                    end
                    
                    figure(ci)
                    plot(dF,'-b')
                    set(gcf,'position',[200 500 1500 400])
                    
                    %%
                    fprintf('Deconvolving ...\n');
                    g = 0.95;
                    [c, s, b, g, lam, active_set] = sc_constrained_oasisAR1(double(dF), g, [], [], [], [], []);
                    % [c, s, b, g, lam, active_set] = constrained_oasisAR1(double(dF), g, [], true, [], [], []);
                    
                    % zero padding after the last spike
                    sp = zeros(1,length(dF));
                    sp(1:length(s)) = s;
                    
                    dFsp.slice(si).cell(ci).dF = dF;
                    dFsp.slice(si).cell(ci).c  = c';
                    dFsp.slice(si).cell(ci).sp = sp;
                    dFsp.slice(si).cell(ci).b  = b;
                    dFsp.slice(si).cell(ci).g  = g;
                    dFsp.slice(si).cell(ci).lam = lam;
                    dFsp.slice(si).cell(ci).neuropilCoefficient = dat.stat(ci).neuropilCoefficient;
                    telapsed = toc(tstart);
                    fprintf('done in %.1f sec.\n',telapsed);
                    
                    S = whos;
                    mem = sum([S.bytes])*1e-6;
                    fprintf('%.1f MB used\n',mem);
                    
                else
                    dFsp.slice(si).cell(ci).sp = [];
                end
                dFsp.slice(si).cell(ci).roi = ci;
            end
        end
    end
    save(fullfile(folder_name,data_file.name),'dFsp','-append')
end
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
        folder_name = fullfile('\\research.files.med.harvard.edu\Neurobio\HarveyLab\Tier1\Shin\ShinDataAll\Suite2P',mouseID,num2str(date_num));
        try
            ops_file = sprintf('regops_%s_%d',mouseID,date_num);
            temp = load(fullfile(folder_name,ops_file));
        catch
            [ops_file,folder_name] = uigetfile([folder_name,'regops*.mat'],'MultiSelect','off');
            temp = load(fullfile(folder_name,ops_file));
        end
        ops0 = temp.ops1{1};
        % data_file = sprintf('F_%s_%d_plane1_proc.mat',mouseID,date_num);
        data_file = dir(fullfile(folder_name,'*_proc.mat'));
    end
    
    try
        movFile = fullfile(ops0.RootDirRaw,'FOV1_001_001.tif'); % SI4?
    catch
        movFile = fullfile(ops0.RootDirRaw,'FOV1_00001_00001.tif'); % SI2016
    end
    
    if isempty(gcp('nocreate'))
        parpool(4);
    end
    
    tstart = tic;
    dat = [];
    msg = [];
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
            num_roi = length(dat.stat);
            
            fit_mode = 'prctile'; %'exp_linear';
            base_prctile = 50;
            
            if base_prctile==50
                fprintf('Using the median F for baseline.\n\n')
            end
            
            %% copying F file with a new name
            file_name_from = sprintf('F_%s_%d_plane1_proc.mat',mouseID,date_num);
            file_name_to   = sprintf('F_%s_%d_plane1_proc_%dprctile.mat',mouseID,date_num,base_prctile);
            [success,~,~] = copyfile(fullfile(folder_name,file_name_from),fullfile(folder_name,file_name_to));

            if success
                fprintf('%s was copied to %s successfully.\n\n',file_name_from,file_name_to);
            else
                error('%s was not copied appropriately.\n\n',file_name_from);
            end
            %%            
            fprintf('Deconvolving F\n\n')
            tstart = tic;
            ROI_num = find([dat.stat(1:num_roi).iscell]);
            
            % check if neuropil subtraction was successful
            if any(isnan([dat.stat.neuropilCoefficient]))
                if all(isnan([dat.stat.neuropilCoefficient]))
                    msg = addMsg(msg,sprintf('neuropil subtraction failed in all cells.'));
                else
                    msg = addMsg(msg,sprintf('neuropil subtraction failed in some cells.'));
                end
            end
            
            for ri = 1:num_roi
                if isgood(ri)
                    %%
                    fprintf('Analyzing ROI %d (%d / %d cells)\n',ri,sum(isgood(1:ri)),sum(isgood));
                    fprintf('Computing dF/F ...... ')
                    if 0
                        % compute dF
                        dF = Fcell(ri,:) - bsxfun(@times, FcellNeu(ri,:), dat.stat(ri).neuropilCoefficient); % using the method in extractSignals.m
                    else
                        % compute dF/F
                        rawTrace = Fcell(ri,:);
                        if ~isnan(dat.stat(ri).neuropilCoefficient)
                            subTrace = Fcell(ri,:) - bsxfun(@times, FcellNeu(ri,:), dat.stat(ri).neuropilCoefficient); % using the method in extractSignals.m
                        else
                            % if neuropilCoefficient was not properly
                            % computed, use rawTrace;
                            subTrace = rawTrace;
                        end
                        rawTraceBase = getF_(rawTrace,fit_mode,ops0.imageRate,base_prctile);
                        subTraceBase = getF_(subTrace,fit_mode,ops0.imageRate,base_prctile);
                        dF = (subTrace - subTraceBase)./rawTraceBase;
                    end
                    fprintf('Done\n')
                    
                    % figure(ri)
                    % plot(dF,'-b')
                    % set(gcf,'position',[200 500 1500 400])
                    
                    %%
                    fprintf('Deconvolving ...... ');
                    g = 0.95;
                    [c, s, b, g, lam, ~] = sc_constrained_oasisAR1(double(dF), g, [], [], [], [], []);
                    % [c, s, b, g, lam, active_set] = constrained_oasisAR1(double(dF), g, [], true, [], [], []);
                    
                    % zero padding after the last spike
                    sp = zeros(1,length(dF));
                    sp(1:length(s)) = s;
                    
                    DF{ri} = dF;
                    C{ri} = c';
                    SP{ri} = sp;
                    B{ri} = b;
                    G{ri} = g;
                    LAM{ri} = lam;
                    NC{ri} = dat.stat(ri).neuropilCoefficient;
                    
                    % dFsp.slice(si).roi(ri).dF = dF;
                    % dFsp.slice(si).roi(ri).c  = c';
                    % dFsp.slice(si).roi(ri).sp = sp;
                    % dFsp.slice(si).roi(ri).b  = b;
                    % dFsp.slice(si).roi(ri).g  = g;
                    % dFsp.slice(si).roi(ri).lam = lam;
                    % dFsp.slice(si).roi(ri).neuropilCoefficient = dat.stat(ri).neuropilCoefficient;
                    % telapsed = toc(tstart);
                    % fprintf('Done in %.1f sec.\n',telapsed);
                    
                    %S = whos;
                    %mem = sum([S.bytes])*1e-6;
                    %fprintf('%.1f MB used\n\n',mem);
                else
                    % dFsp.slice(si).roi(ri).sp = [];
                    % SP{ri} = [];
                end
                % dFsp.slice(si).roi(ri).roi = ri;
                % telapsed = toc(tstart);
                % fprintf('%d sec elapsed\n',si,round(telapsed))
            end
        end
    end
    
    % store the computed variables in a structure
    for si = 1
        for ri = 1:num_roi
            if isgood(ri)
                dFsp.slice(si).roi(ri).dF = DF{ri};
                dFsp.slice(si).roi(ri).c  = C{ri};
                dFsp.slice(si).roi(ri).sp = SP{ri};
                dFsp.slice(si).roi(ri).b  = B{ri};
                dFsp.slice(si).roi(ri).g  = G{ri};
                dFsp.slice(si).roi(ri).lam= LAM{ri};
                dFsp.slice(si).roi(ri).neuropilCoefficient = NC{ri};
                dFsp.slice(si).roi(ri).roi = ri;
            end
        end
    end
    
    fprintf('Saving %s ... ',file_name_to);
    save(fullfile(folder_name,file_name_to),'dFsp','-append')
    fprintf('Done\n\n');
    
    %% save Lite file
    Nframes = dat.ops.Nframes;
    temp = nan(length(ROI_num),Nframes);
    
    for ci = 1:length(ROI_num)
        temp(ci,:) = dFsp.slice.roi(ROI_num(ci)).sp;
    end
    
    clear dFsp
    dFsp = temp; % dFsp just contains spike data
    [~,name,ext] = fileparts(file_name_to);
    file_name_to_lite = [name,'_Lite',ext];
    fprintf('Saving %s ... ',file_name_to_lite);
    save(fullfile(folder_name,file_name_to_lite),'dFsp','ROI_num','msg');
    fprintf('Done\n\n');
    
end

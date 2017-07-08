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
        % data_file = sprintf('F_%s_%d_plane1_proc.mat',mouseID,date_num);
        data_file = dir(fullfile(folder_name,'*_proc.mat'));
    end
    
    movFile = fullfile(ops0.RootDirRaw,'FOV1_001_001.tif'); % SI4?
    if ~exist(movFile,'file')
        movFile = fullfile(ops0.RootDirRaw,'FOV1_00001_00001.tif'); % SI2016
    end
    
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
            for ci = startCellInd:size(Fcell,1)
                if isgood(ci)
                    %%
                    fprintf('Analyzing ROI %d (%d / %d cells)\n',ci,sum(isgood(1:ci)),sum(isgood));
                    fprintf('Computing dF/F ...... ')
                    if 0
                        % compute dF
                        dF = Fcell(ci,:) - bsxfun(@times, FcellNeu(ci,:), dat.stat(ci).neuropilCoefficient); % using the method in extractSignals.m
                    else
                        % compute dF/F
                        rawTrace = Fcell(ci,:);
                        subTrace = Fcell(ci,:) - bsxfun(@times, FcellNeu(ci,:), dat.stat(ci).neuropilCoefficient); % using the method in extractSignals.m
                        rawTraceBase = getF_(rawTrace,fit_mode,ops0.imageRate,base_prctile);
                        subTraceBase = getF_(subTrace,fit_mode,ops0.imageRate,base_prctile);
                        dF = (subTrace - subTraceBase)./rawTraceBase;
                    end
                    fprintf('Done\n')
                    %{
                    figure(ci)
                    plot(dF,'-b')
                    set(gcf,'position',[200 500 1500 400])
                    %}
                    
                    %%
                    fprintf('Deconvolving ...... ');
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
                    fprintf('Done in %.1f sec.\n',telapsed);
                    
                    S = whos;
                    mem = sum([S.bytes])*1e-6;
                    fprintf('%.1f MB used\n\n',mem);
                else
                    dFsp.slice(si).cell(ci).sp = [];
                end
                dFsp.slice(si).cell(ci).roi = ci;
            end
        end
    end
    save(fullfile(folder_name,file_name_to),'dFsp','-append')
end
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
                    
                    fprintf('Analyzing ROI %d (%d / %d cells)\n',ci,sum(isgood(1:ci)),sum(isgood));
                    fprintf('Computing Neuropil coefficient ... ')
                    % compute Neuropil coefficient (as in Acq2P)from plotNeuropilTraces.m
                    fBody = Fcell(ci,:);
                    fNeuropil = FcellNeu(ci,:);

                    % Remove excluded frames (removing them seems to be the most
                    % acceptable solution, since many functions below don't deal well
                    % with nans, and interpolation might skew results):
                    % fBody(excludeFrames) = [];
                    % fNeuropil(excludeFrames) = [];

                    % Remove bleaching:
                    % fBody = deBleach(fBody, 'runningAvg',9001);
                    % fNeuropil = deBleach(fNeuropil, 'runningAvg',9001);
                    [fBody, baselineStats] = deBleach(fBody, 'custom_wfun');
                    fNeuropil = deBleach(fNeuropil, 'custom_wfun');
                    f0Body = prctile(fBody,10);
                    f0Neuropil = prctile(fNeuropil,10);

                    % Smooth traces:
                    smoothWin = gausswin(smoothWindow)/sum(gausswin(smoothWindow));
                    fBody = conv(fBody, smoothWin, 'valid');
                    fNeuropil = conv(fNeuropil, smoothWin, 'valid');

                    % Extract subtractive coefficient btw cell + neuropil and plot
                    % cutoffFreq = 1; %Cutoff frequency in seconds
                    % a = framePeriod / cutoffFreq;
                    % fBodyHighpass = filtfilt([1-a a-1],[1 a-1], fBody);
                    % fNeuropilHighpass = filtfilt([1-a a-1],[1 a-1], fNeuropil);
                    fBodyHighpass = fBody-f0Body;
                    fNeuropilHighpass = fNeuropil-f0Neuropil;
                    df = smooth(abs(diff(fBodyHighpass)), round(2/framePeriod));
                    isFChanging = df>2*mode(round(df*100)/100);

                    traceSubSelection = ~isFChanging;

                    %nSmooth = numel(smoothWin);
                    %traceSubSelection = baselineStats.w(floor(nSmooth/2):end-1-(nSmooth-floor(nSmooth/2)))==1;

                    neuropilCoef = robustfit(fNeuropilHighpass(traceSubSelection),...
                        fBodyHighpass(traceSubSelection),...
                        'bisquare',2);
                    
                    fprintf('neuropilCoef = %.3f\n',neuropilCoef(2));
                    % from extractROIsBin.m
                    traceCell = Fcell(ci,:);
                    rawTrace = traceCell;
                    traceNeuropil = FcellNeu(ci,:);
                    subTrace = traceCell - traceNeuropil*neuropilCoef(2); % subCoef = neuropilCoef(2)
                    
                    % from dFcalc.m
                    if 0
                        dF = (subTrace - getF_(subTrace,'custom_wfun'))...
                            ./getF_(rawTrace,'custom_wfun');
                    else
                        dF = subTrace;
                    end
                    
                    if 0
                        figure(1);clf
                        subplot(4,1,1);plot(rawTrace);title('F cell');
                        subplot(4,1,2);plot(traceNeuropil);title('F neuropil');
                        subplot(4,1,3);plot(subTrace);title('F neuropil subtracted');
                        subplot(4,1,4);plot(dF);title('dF');
                        set(gcf,'position',[200 200 1500 1000])
                        drawnow;
                    end
                                        
                    % convert to dF/F without neuropil subtraction
                    % baseFf = getF_(dF(ci,:), 'custom_wfun');
                    % dF(ci,:) = dF(ci,:)./baseFf;
                    % dFneu(ci,:) = dFneu(ci,:)./baseFf;
                    
                    fprintf('Deconvolving ...\n');
                    g = 0.95;
                    [c, s, b, g, lam, active_set] = constrained_oasisAR1(double(dF), g, [], [], [], [], []);
                    
                    dFsp.slice(si).cell(ci).dF = dF;
                    dFsp.slice(si).cell(ci).c  = c';
                    dFsp.slice(si).cell(ci).sp = s';
                    dFsp.slice(si).cell(ci).b  = b;
                    dFsp.slice(si).cell(ci).g  = g;
                    dFsp.slice(si).cell(ci).lam = lam;
                    dFsp.slice(si).cell(ci).neuropilCoef = neuropilCoef;
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
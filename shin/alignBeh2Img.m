function alignBeh2Img(mouse_num,date_num,varargin)

% Align Virmen behavioral data to ScanImage sync pulses

varargin2V(varargin);
mouseID = getMouseID(mouse_num);
if ~exist('ops_file','var')
    folder_name = fullfile('\\research.files.med.harvard.edu\Neurobio\HarveyLab\Tier1\Shin\ShinDataAll\Suite2P',mouseID,num2str(date_num));
    try 
        ops_file = dir(fullfile(folder_name,'regops*.mat'));
        if length(ops_file)>1
            error('Two or more regops files exist. Choose the appropriate file.')
        end
        load(fullfile(folder_name,ops_file.name));
        ops = ops1{1};
    catch
        warning('Could not find the appropriate ops_file.')
    end
    
else
    [folder_name,ops_file] = fileparts(ops_file);
    load(fullfile(folder_name,ops_file));
    ops = ops1{1};
end

nSlices = ops.nplanes;
if nSlices==1
    fastZDiscardFlybackFrames = 0;
    sliceSet = 1;
else
    fastZDiscardFlybackFrames = 1;
end
nSlices = nSlices - fastZDiscardFlybackFrames;

if ~exist('acqName','var')
    acqName = 'FOV1';
end

if ~exist('sliceNum','var')
    sliceNum = 1;
end

if ~exist('channelNum','var')
    channelNum = 1;
end

chunk_size = 1e6;

temp = userpath;
matlab_dir = temp(1:end-1);
mat_file_name = sprintf('%s%s\\%s_%s',...
    '\\research.files.med.harvard.edu\Neurobio\HarveyLab\Tier1\Shin\ShinDataAll\Current Mice\',...
    ops.mouse_name,ops.mouse_name,ops.date);
mat_file_name = [mat_file_name,'.mat'];
load(mat_file_name);
m_time_stamp = data(1,:);

switch vr.computerID
    case 3
        MatlabPulseMode = 'analog';
    case 4
        MatlabPulseMode = 'digital';
end

nFramesTotal = ops.Nframes * ops.nplanes; % default num of frames

switch MatlabPulseMode
    case 'analog' % recorded by Wavesurfer
        ImgPath = sprintf('\\\\research.files.med.harvard.edu\\Neurobio\\HarveyLab\\Tier2\\Shin\\ShinDataAll\\Imaging\\%s\\%s\\',ops.mouse_name,ops.date);
        obj_file = sprintf('%s_%s_FOV1_00001.mat',ops.mouse_name,ops.date);
        load(fullfile(ImgPath,obj_file));
        eval(['obj = ',obj_file(1:end-4),';']);
        samp_rate = 1e3;
        file_name = [obj.defaultDir,'FOV1_0001.h5'];
        file_name = changePath4Server(file_name);
        file_name = strrep(file_name,'HarveyLab\Shin','HarveyLab\Tier2\Shin');
        try
            wsData = h5read(file_name,'/sweep_0001/analogScans');
        catch
            wsData = h5read(file_name,'/sweep_0002/analogScans');
        end

        SIsig = wsData(:,4)';
        temp = SIsig>7800;
        SI_sig_rise = find(diff(temp)==1)+1;
        SI_sig_fall = find(diff(temp)==-1)+1;

        NIsig = double(wsData(:,5))'*0.15/462.7; % matlab pulse -- convert to voltage
        N = length(NIsig);
        num_chunk = ceil(N/chunk_size);

        n = length(m_time_stamp);
        sig_rise = nan(1,n);
        sig_fall = nan(1,n);
        sig_rise_t = nan(1,n);
        amp = nan(1,n);

        sig_fall_temp = 0;
        split_flag = 0;
        chunk_end = 0;

        x = [[0,0,NIsig(1,1:end-2)];NIsig(1,:)];
        dx = diff(x,1,1);
        threshold = 0.07;
        rise_all = find(dx>threshold);
        fall_all = find(dx<-threshold);
        m = 1;
        y = nan(1,1e5);
        j = 1;
        k = 1;
        for i = 1:num_chunk
            chunk_end(i+1) = chunk_end(i) + chunk_size;
            a = rise_all(find(rise_all<=chunk_end(i+1)+1,1,'last'));
            b = fall_all(find(fall_all<=chunk_end(i+1),1,'last'));
            while isempty(a) || a>b || b-a<2
                chunk_end(i+1) = chunk_end(i+1) + chunk_size;
                a = rise_all(find(rise_all<=chunk_end(i+1)+1,1,'last'));
                b = fall_all(find(fall_all<=chunk_end(i+1),1,'last'));
            end

            if chunk_end(i+1)>N
                sig = NIsig(chunk_end(i)+1:end);
            else
                sig = NIsig(chunk_end(i)+1:chunk_end(i+1));
            end

            x = [[0,0,sig(1:end-2)];sig];
            dx = diff(x,1,1);
            rise = find(dx>threshold);
            fall = find(dx<-threshold);
            sig_fall_temp = 0;

            while 1

                sig_rise_temp = rise(find(rise>sig_fall_temp,1,'first'));

                if k==1
                    sig_rise_temp = rise(find(rise>sig_rise_temp+1,1,'first'));
                end
                if isempty(sig_rise_temp)
                    break
                end
                sig_fall_temp = fall(find(fall>sig_rise_temp,1,'first'));
                if isempty(sig_fall_temp)
                    error('sig_fall_temp should not be empty!!')
                end

                sig_rise_temp2 = rise(find(rise>sig_fall_temp,1,'first'));
                if isempty(sig_rise_temp2)
                    sig_rise_temp2 = sig_fall_temp + 10;
                end
                
                if k>1 % skip k=59950 (VS045 170911) k=105652 (VS045 170913) k=28779 (VS045 170922) k=9852 (VS035 170508)
                    if sig_fall_temp - sig_rise_temp < 2
                             % transient artifact (spike) at the baseline
                        if  abs(sig(sig_fall_temp+2) - sig(sig_rise_temp-2)) < 0.01
                            continue
                        else % transient artifact (spike) at the pulse onset
                            sig_fall_temp = fall(find(fall>sig_fall_temp+1,1,'first'));
                        end
                    end
                end
                    % transient artifact (dip) on the pulse. Look for the next fall.
                    % Be aware!! there is a slight chance that the dip could 
                    % occurs at the end of chunks.
                if sig_rise_temp2 - sig_fall_temp < 2
                    sig_fall_temp = fall(find(fall>sig_fall_temp+1,1,'first'));
                end

                sig_rise(k) = sig_rise_temp + chunk_end(i);
                sig_fall(k) = sig_fall_temp + chunk_end(i);
                sig_rise_t(k) = (sig_rise(k) - sig_rise(1))/30e3;
                amp(k) = max(sig(sig_rise_temp:sig_fall_temp));

                if mod(k,1e4)==0
                    fprintf('%d pulses detected\n',k)
                end
                j = j+1;
                k = k+1;
            end
        end
        NI_time_stamp = sig_rise';
        % NI_time_stamp = (sig_rise' - sig_rise(1));
        num_pulse = find(isfinite(amp),1,'last');
        
        if num_pulse == length(m_time_stamp)
            fprintf('All Virmen pulses were detected!\n')
        elseif num_pulse < length(m_time_stamp)
            warning('Some of Virmen sync pulses were not detected');
        end
        
        data = data(:,1:num_pulse);
        
    case 'digital' % recorded by Matthias' Sync obj
        % detecting MATLAB sync pulses
        file_name = sprintf('%s%s\\%s\\%s_0001.h5','\\research.files.med.harvard.edu\Neurobio\HarveyLab\Tier2\Shin\ShinDataAll\Imaging\',ops.mouse_name,ops.date,'FOV1');
        if exist(file_name,'file')
            a = ws.loadDataFile(file_name);
            samp_rate = 2e3;
            for i = 1:length(a.header.Acquisition.ActiveChannelNames)
                pick_NI(i) = ~isempty(strfind(a.header.Acquisition.ActiveChannelNames{i},'VirmenSync'));
                pick_SI(i) = ~isempty(strfind(a.header.Acquisition.ActiveChannelNames{i},'ScanImageSync'));
            end
            % detecting Virmen sync pulses
            NIsig = a.sweep_0001.analogScans(:,pick_NI);
            % detecting ScanImage sync pulses
            SIsig = a.sweep_0001.analogScans(:,pick_SI);
        end

        file_name = dir(sprintf('%s%s\\%s\\%s_syncData*.mat','\\research.files.med.harvard.edu\Neurobio\HarveyLab\Tier2\Shin\ShinDataAll\Imaging\',ops.mouse_name,ops.date,'FOV1'));
        if ~isempty(file_name)
            file_name = sprintf('%s%s\\%s\\%s','\\research.files.med.harvard.edu\Neurobio\HarveyLab\Tier2\Shin\ShinDataAll\Imaging\',ops.mouse_name,ops.date,file_name.name);
            if exist(file_name,'file')
                temp = load(file_name);
                if isfield(temp,'sync')
                    s = temp.sync;
                else
                    s = temp.s;
                end
                num_channel = length(s.chNames);
                samp_rate = s.prop.snDataLogProperties.Rate;
                for i = 1:num_channel
                    pick_NI(i) = ~isempty(strfind(s.chNames{i},'VirmenSync'));
                    pick_SI(i) = ~isempty(strfind(s.chNames{i},'ScanImageSync'));
                end
                % detecting Virmen sync pulses
                NIsig = s.ch(:,pick_NI);
                % detecting ScanImage sync pulses
                SIsig = s.ch(:,pick_SI);
            end
        end
        temp = NIsig>2.5;
        NI_sig_rise = find(diff(temp)==1)+1;
        NI_sig_fall = find(diff(temp)==-1)+1;
        
        temp = SIsig>2.5;
        SI_sig_rise = find(diff(temp)==1)+1;
        SI_sig_fall = find(diff(temp)==-1)+1;
        % si_time_stamp = find(diff(SIsig)>2.6)*1e3/samp_rate; % scan image pulse
        
        if NI_sig_rise(1) < NI_sig_fall(1)
            rise_first = true;
            NI_time_stamp = sort([NI_sig_rise;NI_sig_fall]);
        else
            rise_first = false;
            % This indicates (i) successful reset OR (ii) failed reset 
            if length(NI_sig_rise)+length(NI_sig_fall) == length(m_time_stamp)+1
                % case (i): successful reset, the first sig_fall occurred before
                % the 1st Virmen iteration.
                NI_time_stamp = sort([NI_sig_rise;NI_sig_fall]);
                NI_time_stamp = NI_time_stamp(2:end);

            elseif length(NI_sig_rise)+length(NI_sig_fall) == length(m_time_stamp)-1
                % case (ii): failed reset, the first sig_fall occurred at the end
                % of the 2nd Virmen iteration.
                NI_time_stamp = sort([NI_sig_rise;NI_sig_fall]);
                data = data(:,2:end);
            end
        end

        if (length(NI_sig_rise) + length(NI_sig_fall)) < length(m_time_stamp)
            warning('Some of Virmen sync pulses were not detected');
        end
end

if length(SI_sig_rise) > nFramesTotal * (nSlices + fastZDiscardFlybackFrames)
    warning('The number of ScanImage sync pulses exceeds the saved number of image frames');
    init_ind  = find(diff(SI_sig_rise)>40e-3*samp_rate)+1;
    init_ind = [1,init_ind,length(SI_sig_rise)+1];
    k = 1;
    minBlockSize = 1e4;
    for i = 1:length(init_ind)-1
        pick = init_ind(i):init_ind(i+1)-1;
        if length(init_ind(i):init_ind(i+1)-1) >= minBlockSize;
            blockFrames{k} = pick;
            k = k+1;
        end
    end
    block_ind = uiSelectBlock(SIsig(1:1e3:end),blockFrames);
    pick = blockFrames{block_ind};
    SI_sig_rise = SI_sig_rise(pick);
    SI_sig_rise = SI_sig_rise(1:nFramesTotal);
end

if iscolumn(NI_time_stamp)
    NI_time_stamp = NI_time_stamp';
end
if iscolumn(SI_sig_rise)
    SI_sig_rise = SI_sig_rise';
end

% convert time-stamp indices to ms
NI_time_stamp_ms = NI_time_stamp*1e3/samp_rate;
SI_time_stamp_ms = SI_sig_rise*1e3/samp_rate;

for i = 1:size(data,1)
    temp = interp1(NI_time_stamp_ms,data(i,:),SI_time_stamp_ms);
    interpData(i,:) = temp;
end

BM = [SI_time_stamp_ms;interpData];

save_name = sprintf('B_%s_%s.mat',ops.mouse_name,ops.date);
save(fullfile(folder_name,save_name),'BM');

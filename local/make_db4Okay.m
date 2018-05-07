
github_dir = 'C:\Users\Shin\Documents\GitHub\DmtsShared\DMTS';
load(fullfile(github_dir,'info_set_good_okay.mat'));

Suite2P_dir = '\\research.files.med.harvard.edu\Neurobio\HarveyLab\Tier1\Shin\ShinDataAll\Suite2P';
Imaging_dir = '\\research.files.med.harvard.edu\Neurobio\HarveyLab\Tier2\Shin\ShinDataAll\Imaging';

load(fullfile(github_dir,'MotionCorrectInfo.mat'));

mouse_set = [9,15,20,31,35,45];
for mi = 1:length(mouse_set)
    pick = cell2mat(info_set(:,1))==mouse_set(mi);
    date_set{mi} = cell2mat(info_set(pick,2));
    info_set_mouse{mi} = info_set(pick,:);
end

i = 0;
for mi = 1:length(mouse_set)
    
    if isempty(date_set{mi})
        continue
    end
    
    mouse_num = mouse_set(mi);
    mouseID = getMouseID(mouse_num);

    
    for di = 1:length(date_set{mi})
        
        date_num = date_set{mi}(di);        
        ind = find([info_set_mouse{mi}{:,2}]==date_num);
        FOV_info = info_set_mouse{mi}(ind,:);
        F_file_name = sprintf('F_%s_%d_plane1.mat',mouseID,date_num);
        
        skip_flag = false;
        
        % skip if F file already exists
        F_file_path = fullfile(Suite2P_dir,mouseID,num2str(date_set{mi}(di)),F_file_name);
        file_info = dir(F_file_path);
        if exist(F_file_path,'file')
            skip_flag = true;
        end
        
        % but do not skip if old motion correction was used to create F file
        pick = cell2mat(MotionCorrectInfo(:,1))==mouse_num & cell2mat(MotionCorrectInfo(:,2))==date_num;
        if ~strcmp(MotionCorrectInfo{pick,3},'lucasKanade_plus_nonrigid')
            skip_flag = false;
        end
        
        % skip if img files have not been motion corrected
        if 0 %isempty(MotionCorrectInfo{di,3})
            skip_flag = true;
        end
        
        % skip if F file was create on 18/04/26 or later
        % if datenum(file_info.date) > datenum('2018/04/26')
        %     skip_flag = true;
        % end
        
        if FOV_info{2}==160607
            skip_flag = true;
        end

        if skip_flag
            % fprintf('skipping %s %d\n',mouseID,date_num);
            continue
        end        
        
        i = i+1;
                
        db(i).mouse_name    = mouseID;
        db(i).date          = num2str(FOV_info{2});
        db(i).area          = FOV_info{3};
        db(i).depth         = num2str(FOV_info{4});
        db(i).expts         = [1];
        db(i).nchannels     = 1;
        db(i).gchannel      = 1;
        db(i).nplanes       = 1;
        db(i).diameter      = 16; % 12 is too small, 20 is too large
        db(i).comments      = sprintf('%s %s %s %sum',db(i).mouse_name,db(i).date,db(i).area,db(i).depth);
        
        fprintf('%s %s %s %sum\n',db(i).mouse_name,db(i).date,db(i).area,db(i).depth);
    end    
end
db = db(end-2);

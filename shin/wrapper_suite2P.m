function wrapper_suite2P

    mouse_set = [20,15,9];
    date_set = {[161114,161117,161120,161221],[161129,161227,170112,161205,161221],[160525,160516,160608,160613]};
    
    for mi = 1:length(mouse_set)
        for di = 1:length(date_set{mi})
            mouse_num = mouse_set(mi);
            initials = getInitials(mouse_num);
            mouseID = sprintf('%s%03d',initials,mouse_num);
            date_num = num2str(date_set{mi}(di));
            fprintf('Processing %s %s',mouseID,date_num)
            
            folder_name = '\\research.files.med.harvard.edu\Neurobio\HarveyLab\Shin\ShinDataAll\Suite2P\';
            ops_file = sprintf('regops_%s_%s.mat',mouseID,date_num);
            full_ops_file = fullfile(folder_name,mouseID,date_num,ops_file);
            save_file = sprintf('B_%s_%s.mat',mouseID,date_num);
            full_save_file = fullfile(folder_name,mouseID,date_num,save_file);
            if ~exist(full_ops_file,'file')
                fprintf('... skipping:\tNo Suite2P ops_file.\n')
                continue
            elseif exist(full_save_file,'file')
                fprintf('... skipping:\tB_data already exists.\n')
                continue
            end
            alignBeh2Img('ops_file',full_ops_file);
            fprintf('... done.\n')
        end
    end

return
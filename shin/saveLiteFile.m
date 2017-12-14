function saveLiteFile(mouse_num,date_num)
    
    mouseID = getMouseID(mouse_num);
    dir_from = fullfile('\\research.files.med.harvard.edu\Neurobio\HarveyLab\Tier1\shin\ShinDataAll\Suite2P',mouseID,num2str(date_num));
    dir_to   = fullfile('\\research.files.med.harvard.edu\Neurobio\HarveyLab\Tier1\shin\ShinDataAll\Suite2P',mouseID,num2str(date_num));
    
    if ~exist(dir_to,'dir')
        mkdir(dir_to);
        fprintf('Created a new directory:\n%s\n',dir_to)
    end

    file_name_from = sprintf('F_%s_%d_plane1_proc_50prctile.mat',mouseID,date_num);
    [~,file_name,file_ext]   = fileparts(file_name_from);
    file_name_to   = [file_name,'_Lite',file_ext];
    
    fprintf('%Saving %s\n\n',file_name_to);
    
    load(fullfile(dir_from,file_name_from))
    
    cell_set = find([dat.stat.iscell]);
    Nframes = dat.ops.Nframes;
    temp = nan(length(cell_set),Nframes);
    
    for ci = 1:length(cell_set)
        temp(ci,:) = dFsp.slice.cell(cell_set(ci)).sp;
        ROI_num(ci) = cell_set(ci);
    end
    
    clear dFsp
    dFsp = temp; % dFsp just contains spike data
    
    save(fullfile(dir_to,file_name_to),'dFsp','ROI_num');
 
end
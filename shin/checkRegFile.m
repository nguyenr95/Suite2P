function checkRegFile
    
    nimgbatch = 2000;
    ix = 0;

    fclose all;
    ResultsSavePath = 'C:\Users\Shin\Documents\MATLAB\ShinDataAll\Suite2P\LT009\160520\';
    load(fullfile(ResultsSavePath,'regops_LT009_160520_plane1.mat'))
    ops.RegFileRoot = 'E:\Imaging\Suite2P\';
    ops.RegFile = fullfile(ops.RegFileRoot,'tempreg_plane1.bin');
    [LyU, LxU] = size(ops.mimg);
    
    fid = fopen(ops.RegFile, 'r');

    tic

    % while 1
        mov = fread(fid,  LyU*LxU*nimgbatch, '*int16');
        if isempty(mov)
            % break;
        end
        mov = reshape(mov, LyU, LxU, []);
        mov = mov(ops.yrange, ops.xrange, :);
        mov = single(mov);
        NT= size(mov,3);

        % mov = single(reshape(mov, [], NT));
        % fprintf('Frame %d done in time %2.2f \n', ix, toc)
    % end
    fclose(fid);
    tiffWrite(mov,'test',ops.RegFileRoot, 'int16');
end
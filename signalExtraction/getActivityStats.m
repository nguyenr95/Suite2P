function [stat, F, Fneu] = getActivityStats(ops, stat, F, Fneu)

indNoNaN    = find(~ops.badframes);
ix          = cumsum(~ops.badframes) + 1;
ix          = ix(ops.badframes);
ix(ix>numel(indNoNaN))  = numel(indNoNaN);

F(:, ops.badframes)  = F(:,    indNoNaN(ix));
Fneu(:, ops.badframes)  = Fneu(:, indNoNaN(ix));

ops.fs           = getOr(ops, 'fs', ops.imageRate/ops.nplanes);
%
try 
    % figure out the ICA coefficients here <--- FIX THIS
    % my_ica was commented out in refs/remotes/cortex-lab/master on
    % 17/05/17, keep updated with the associated changes!
    [coefNeu, inomax]   = my_ica(F', Fneu', ops.fs, 0.7);

    coefNeu = 0.7 * ones(1, size(F,1));

    dF                  = F - bsxfun(@times, Fneu, coefNeu(:));

    % dF          = F - Fneu;

    sd           = std(dF, [], 2);
    sdN          = std(Fneu, [], 2);

    sk(:, 1) = skewness(dF, [], 2);
    sk(:, 2) = sd./sdN; 
    sk(:, 3) = (max(dF, [], 2)-median(dF, 2))./sd;
    sk(:, 4) = (prctile(dF, 95, 2)-median(dF, 2))./sd;

    for j = 1:numel(stat)
        stat(j).dFstat           = sk(j,:);
        stat(j).skew             = sk(j,1);
        stat(j).std              = sk(j,2);
        stat(j).maxMinusMed      = sk(j,3);
        stat(j).top5pcMinusMed   = sk(j,4);
        stat(j).blockstarts      = [0 cumsum(ops.Nframes)];
        stat(j).iplane                 = ops.iplane;
        stat(j).neuropilCoefficient    = coefNeu(j); 
    end
catch
    save(sprintf('%s/FandFneu_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, ops.iplane), 'F', 'Fneu', '-v7.3')
    for j = 1:numel(stat)
        stat(j).dFstat           = NaN;
        stat(j).skew             = NaN;
        stat(j).std              = NaN;
        stat(j).maxMinusMed      = NaN;
        stat(j).top5pcMinusMed   = NaN;
        stat(j).blockstarts      = [0 cumsum(ops.Nframes)];
        stat(j).iplane                 = ops.iplane;
        stat(j).neuropilCoefficient    = NaN; 
    end
end

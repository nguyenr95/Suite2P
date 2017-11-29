function h = classROI(h)

[~, filename1] = fileparts(h.dat.cl.fpath);
set(h.text50,'String', [filename1 '.mat']);

% st = h.st;
% prior = h.prior;
% statLabels = h.statLabels;

% this load st, prior and statLabels

computerName = getComputerName;

try
    load(h.dat.cl.fpath);
catch
    [pathstr,name,ext] = fileparts(h.dat.cl.fpath);
    switch computerName
        case 'shin-pc'
            pathstr = 'C:\Users\Shin\Documents\GitHub\Suite2P\configFiles\';
            h.dat.cl.fpath = [pathstr,name,ext];
        case 'harveylab41223'
            pathstr = 'C:\Users\harveylab\Documents\GitHub\Suite2P\configFiles\';
            h.dat.cl.fpath = [pathstr,name,ext];
    end
    load(h.dat.cl.fpath);
end
stat = h.dat.stat;

ilbl = false(1, numel(statLabels));
% st0 = repmat( mean(st,1), numel(stat), 1);
st0 = zeros( numel(stat), numel(statLabels));
for j = 1:numel(statLabels)
    if isfield(stat, statLabels{j})
        ilbl(j) = true;
        st0(:,j) = eval(sprintf('[stat.%s]', statLabels{j}));
    end
end
st0(isnan(st0)) = 2;
ilbl(1) = false;

out = classDialogBox;

if out
    % model       = buildaLDA(prior, st, statLabels);
    % [Ypred, ps] = evaluateLDA(model, st0);
    model       = buildHist(prior, st, statLabels);
    [Ypred, ps] = evaluateHist(model, st0);

    for j = 1:length(stat)
        stat(j).iscell = Ypred(j) > h.dat.cl.threshold; 
        stat(j).cellProb  = Ypred(j);
        % overwrite manual selections -- SK 17/02/08
        stat(j).iscell(stat(j).manual>1e-3) = 1;
        stat(j).iscell(stat(j).manual<-1e-3) = 0;
    end
end

h.dat.stat      = stat;
h.st0           = st0;
h.statLabels    = statLabels;

% h.dat.cl.rands = 0.1 + .8 * Ypred;    

% set(h.text51,'String', fname);

% ps

function out = classDialogBox

    % Construct a questdlg with three options
    choice = questdlg('(Re-) apply a classifier?', ...
        'Dessert Menu', ...
        'No','Yes','No');
    % Handle response
    switch choice
        case 'Yes'
            out = true;
        case 'No'
            out = false;
    end

end
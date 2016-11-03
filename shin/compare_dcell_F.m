
nFramesTotal = ops.Nframes;
n_cell = 592;
sp = zeros(n_cell,nFramesTotal);

for ci = 1:n_cell
    sp(ci,cl.dcell{ci}.st) = cl.dcell{ci}.c;
end


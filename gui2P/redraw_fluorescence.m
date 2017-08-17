function redraw_fluorescence(h)
axes(h.axes4)
hold off

ichosen = h.dat.F.ichosen;
F = [];
Fneu = [];
for j = 1:numel(h.dat.Fcell)
    F    = cat(2, F, h.dat.Fcell{j}(ichosen, :));
    Fneu = cat(2, Fneu, h.dat.FcellNeu{j}(ichosen, :));
end

F_smooth = my_conv_local(medfilt1(double(F), 3), 3);
plot(F_smooth)
axis tight
hold on

coefNeu = 1;
if isfield(h.dat.stat, 'neuropilCoefficient')
    coefNeu = h.dat.stat(ichosen).neuropilCoefficient;
end

baseline = 0;
if isfield(h.dat.stat, 'baseline')
    baseline = h.dat.stat(ichosen).baseline;
end

if isfield(h.dat, 'FcellNeu')
    Fneu_smooth = coefNeu * my_conv_local(medfilt1(double(Fneu), 3), 3);
    plot(Fneu_smooth)
    plot(F_smooth - Fneu_smooth,'-r')
end

box off
axis off
set(gca, 'xcolor', 'w')
% plot([0 NT], [0 0], 'k', 'Linewidth', 2)
% axis off

%% To find swim bouts as peaks in the ramp stimulation

clear
close all
%%
path1 = '/home/ljp/Science/Projects/Rolling/Analysis/2018-01/2018-01-16/'
folder1 = 'Run.15'
% 
% path = '/home/ljp/Science/Projects/Rolling/Data_balancoir/2018-01/'
% folder = 'ramp_90_500_4dpf_01'
load ([path1 folder1 '/workspace.mat'])
sep = find(path1 == filesep); %find the separations in the directory made by a /
Exp = [path1(sep(end-2):end) folder1]; % defines the name of the folder or directory to name the

cd (fullfile(path1, folder1))
%%
N = 10 % number of data to average with the moving average method of smothing
Tail2e_smooth = smoothdata(Tail2e,'movmean',N);

Swb = diff(Tail2e_smooth)

figure; hold on
yyaxis left
plot(Time (fi:inc:(ff-1)), Swb(fi:inc:(ff-1)), 'color', 'c')

% plot(Time (fi:inc:ff), Tail2e(fi:inc:ff), 'LineStyle','-')

plot(Time (fi:inc:ff), Tail2e_smooth(fi:inc:ff),'LineStyle','-')

ylim ([-25 25])
ylabel ('Tail angle (deg)')

yyaxis right
plot (Time(fi:inc:ff),Motor(fi:inc:ff), 'LineStyle','-')
ylabel ('Absolute body angle (deg)')

ylim ([-90 90])
xlabel ('Time (s)')
% 
legend('Diff(tail)',['Tail 2e smoothen movmean' num2str(N)],'Motor Angle', 'Location','Southoutside')
title ([Exp])

%%

[peak_amp , peak_locs ] = findpeaks(abs(Swb),'MinPeakHeight',1.5,'MinPeakDistance',20);
% peak_locs = peak_locs (1:50)
% peak_amp = peak_amp(1:50)
fig = figure
hold on
yyaxis left
plot(Time (fi:inc:(ff-1)), Swb(fi:inc:(ff-1)))
ylim ([-15 15])
ylabel ('Tail velocity')
yyaxis right
plot (Time(fi:inc:ff),Motor(fi:inc:ff), 'LineStyle','-')
plot(Time(peak_locs), Motor(peak_locs),'*','Color','g');
ylim ([-90 90])
ylabel ('Motor Angle')
title( [Exp])
xlabel ('Time')
legend ('diff(tail)', 'motor ang', 'swimm bout', 'Location', 'southoutside')
%%

savefig ('Swim_bouts.fig')
saveas (fig, 'Swim_bouts.pdf')

%%

Angle_peaks = Motor(peak_locs)
save Angle_peaks

%%
cd ([path1, folder1])
dlmwrite ('Angle_peaks.txt', Angle_peaks, '\t')
%%

save workspace




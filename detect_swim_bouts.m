%% To find swim bouts as peaks in the ramp stimulation

clear
close all
%%
path = '/home/ljp/Science/Projects/Rolling/2017-12/2017-12-14/'
folder = 'run.1'
load ([path folder '/AngleDif.mat'])
load ([path folder '/data.mat'])
sep = find(path == filesep); %find the separations in the directory made by a /
Exp = [path(sep(end-2):end) folder]; % defines the name of the folder or directory to name the

cd (fullfile(path, folder))


%% if parameters ff fi inc are not in the workspace, run:

v = VideoReader('video.avi');


Tf= v.NumberOfFrames;%total number of frames
fi = 1;%initital frame
inc = 1;%increament in the number of frames, if you want to analyze all the sequence it should be = 1
ff = Tf;% final frame, it corresponds to the total lenght of the file list

%% Calculate the derivative to find peaks
N = 10 % number of data to average with the moving average method of smothing
Tail_smooth = smoothdata(angle_dif,'movmean',N); % smoothens the angle_dif vector by averaging 10 points
Time = Data.TimeStamp
Motor = Data.MotorAngle
Swb = diff(Tail_smooth) % calculates the derivative from

figure; hold on
yyaxis left
plot(Time (fi:inc:(ff-1)), Swb(fi:inc:(ff-1)), 'color', 'c')

% plot(Time (fi:inc:ff), Tail2e(fi:inc:ff), 'LineStyle','-')

plot(Time (fi:inc:ff), Tail_smooth(fi:inc:ff),'LineStyle','-')

ylim ([-25 25])
ylabel ('Tail angle (deg)')
yyaxis right
plot (Time(fi:inc:ff),Motor(fi:inc:ff), 'LineStyle','-')
ylabel ('Absolute body angle (deg)')

ylim ([-90 90])
xlabel ('Time (s)')
% 
legend('Diff(tail)',['Tail smoothen movmean' num2str(N)],'Motor Angle', 'Location','Southoutside')
title ([Exp])

%% Find peaks and plot on time

[peak_amp , peak_locs ] = findpeaks(abs(Swb),'MinPeakHeight',1.5,'MinPeakDistance',20);

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


%% Find peaks and plot on index number
[peak_amp , peak_locs ] = findpeaks(abs(Swb),'MinPeakHeight',1.5,'MinPeakDistance',20);

fig = figure
hold on

yyaxis left
plot(Swb(fi:inc:(ff-1)))
ylim ([-15 15])
ylabel ('Tail velocity')
yyaxis right
plot (Motor(fi:inc:ff), 'LineStyle','-')
plot(peak_locs, Motor(peak_locs),'*','Color','g');
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
cd ([path, folder])
dlmwrite ('Angle_peaks.txt', Angle_peaks, '\t')
%%

save ('swimbouts')




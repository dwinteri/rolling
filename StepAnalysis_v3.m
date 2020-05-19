% Analyse step stimulation data

clear all
close all

%% Open Folder and loading Data

path = '/home/ljp/Science/Projects/Rolling/2018-03//2018-03-12/';
folder1= '/Run.02/';

%load data from record
load([path folder1 'data.mat'], 'Data')


% if the video is analyzed for 2 ellipses, use this:
load([path folder1 'AngleDif.mat'], 'angle_dif')

tmp = find(path == filesep); %find the separations in the directory made by a /
Exp = [path(tmp(end-2):end) folder1]; % defines the name of the folder or directory to name the plot

%% Load parameters
load 
filename = ([path folder1 'parameters']);
delimiter = {''};
formatSpec = '%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
param = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;
%%
figure
subplot (2,1,1)
% plot(Data.TimeStamp, Data.TailAngle)
plot(Data.TimeStamp,angle_dif)
subplot (2,1,2)
plot(Data.TimeStamp, Data.MotorAngle)



%% Finding steps

%To find the steps on the motor pattern is necessary to find a way to get
%the max and min on the signal.
%We can do this using the findpeaks tool,

Ncycles = 1;%number of cycles you are considering (it can be 10 to analyze cycle by cycle, or 1 to have all data)

%Renaming the data vectors
MotorAngle  = Data.MotorAngle(1:end/Ncycles);

%To analyze with angle calculated by Hugo's program
%If the angles are more than 180:
Data.TailAngle = (Data.TailAngle - (Data.TailAngle >180) *180)
% 
TailAngle   = Data.TailAngle(1:end/Ncycles);
%For 2 ellipses:
% TailAngle   = angle_dif(1:end/Ncycles);

TimeStamp  = Data.TimeStamp(1:end/Ncycles);
dMotorAngle = diff(MotorAngle);% obtain the derivative of MotorAngle

%Smoothen the signal

%filter motor angle
NavgMotor = 10;% arbitrary number of data index that you want to average for filter them
BM = 1./ones(1,NavgMotor);%to put a dot before an operator applies this operator element by element of the matrix or vector
% B creates a vector, in this case [1 1 1 1 1 1 1 1 1 1]
dMotorAngle_filt = filtfilt(BM,NavgMotor,dMotorAngle);
%apply filter to the derivative MotorAngle vector

%filter tail angle
NavgTail = 20;
BT = 1./ones(1,NavgTail);
TailAngle_filt = filtfilt(BT,NavgTail,TailAngle);

% %plot Filtered Data, not filtered and motor angle
% figure;
% plot(dMotorAngle_filt)
% hold on
% plot(dMotorAngle)
% plot (MotorAngle)

[peak_amp , peak_locs ] = findpeaks(abs(dMotorAngle_filt),'MinPeakHeight',0.1,'MinPeakDistance',100);
% find peaks on the filtered signal, using the absolut value (abs),
% Minimun peak hight and min peak distance defines which height will be the
% minimum to consider the peak as a peak and how much is the lower distance
% to consider between two peaks to count them as peaks.
peak_locs = [ 1 peak_locs ]; %redefine the vector peak locks to start from the first peak

% Plot tracking results, to show where are the peaks
n = 3; %number of raw panels to show in the figure
findStep_fig = figure;
ax1(1) = subplot(n,1,1); %plot in the same figure n raws and 1 column
plot(MotorAngle);hold on
plot(peak_locs, MotorAngle(peak_locs),'*');%plot the peak of the motor angle as function of the index of the peak loc
plot(TailAngle)
ylabel('Motor/Tail Angle')
% ylim ([-20 20])

legend('Motor Angle','Peaks', 'Tail Angle' )
title(Exp,'interpreter','none') %Set the Interpreter property as 'none' so that the string X_1 is displayed in the figure as typed, without making 1 a subscript of X.
% to not read the underscores as subindex
ax1(2) = subplot(n,1,2);
plot(dMotorAngle); hold on
plot(dMotorAngle_filt);
plot(peak_locs , dMotorAngle_filt(peak_locs),'*')
ylabel ('Motor derivative')
legend('Raw', 'Filtered')

ax1(3) = subplot(n,1,3);
plot(TailAngle)
legend('Tail Angle')
linkaxes(ax1,'x') %link all the axis under the name of ax1, only in the x axis
xlabel('FrameIndex')
ylabel('Tail Angle')
%ylim ([-20 20])
xlim([0 inf])

% Verify if the peaks are right and correspond correctly to the steps
%% Save Step Analysis
save workspace
cd([path folder1])
save('StepAnalysis.mat')
savefig(gcf,'findStep_fig.fig')
saveas (gcf,'findStep_fig.pdf')

%% Plot all the steps of one cycle

%Create 2 cells with the time and index for each data point

%Plot Motor Angle of all cycles overlaped to select the right peak locs
figure;
for i = 1:8 % i is the number of cycle to iterate the motor angle
   plot(TimeStamp(peak_locs((i-1)*60+1):peak_locs((i-1)*60+61)) - TimeStamp(peak_locs((i-1)*60+1)), MotorAngle(peak_locs((i-1)*60+1):peak_locs((i-1)*60+61)));hold on;
      %find the las peak_loc index of the first cycle and then substrac 1
   %    legend_info(i) = [int2str(i)];

    
    
end


%% Load specific cycle data into temporal variables (creating cells)

NSteps = 8;
cycle = 7 %Total number of cycles (without the first)

for cycle = 1:7%to eliminate the first cycle
    %Plot for one cycle, all the steps from 1 to 8 (2 to 16 deg)
    for i = 1 : 2*NSteps
        T_tmp{i}   = TimeStamp(peak_locs((2*i-1)+(cycle-1)*60):peak_locs((2*i+1)+(cycle-1)*60)) - TimeStamp(peak_locs((2*i)+(cycle-1)*60));
        % To have all the data point of the different steps (12) at the same time point  (respect to the first cycle)
        %find the las peak_loc index of the first cycle and then substrac 1
        ind{i} = [ peak_locs((2*i-1)+(cycle-1)*60):peak_locs((2*i+1)+(cycle-1)*60) ];
        % To generate index for each data point
    end
    
    
    
    T = [-22 : 0.01 : 22];
    % To display all steps, positives and negatives, on the same time, is
    % necessary to sepate positives and negatives for each step. This is also
    % done by making cells.
    % And in order to be able to average, they must be all the same length, so
    % it's necessary to interpolate all the vectors.
    for i = 1 : NSteps
        % Positive values for TailAngle
        TailAngle_p{i,cycle} = interp1(T_tmp{(2*i-1)}, TailAngle_filt(ind{(2*i-1)}) , T ); % interpolate the tail angle in T points
        %                            X    ,underlying function F(x)   , query points
        % defines odds as positive values
        %vq = interp1(x,v,xq) returns interpolated values of a 1-D function at specific query points using linear interpolation. Vector x contains the sample points, and v contains the corresponding values, v(x). Vector xq contains the coordinates of the query points.
        %If you have multiple sets of data that are sampled at the same point coordinates, then you can pass v as an array. Each column of array v contains a different set of 1-D sample values.
        
        % Negative values for TailAngle
        TailAngle_n{i,cycle} = interp1(T_tmp{(2*i)}  , TailAngle_filt(ind{(2*i)})   , T );
        % Positives values for MotorAngle
        MotorAngle_p{i,cycle} = interp1(T_tmp{(2*i-1)}, MotorAngle(ind{(2*i-1)}) , T );
        % Negative values for MotorAngle
        MotorAngle_n{i,cycle} = interp1(T_tmp{(2*i)}  , MotorAngle(ind{(2*i)})   , T );
    end
   
    %Plot in a differente window every cycle 
    
    c = jet(NSteps); % color jet scale for the different steps
    
    Step_fig = figure; %create a figure named Step_fig
    ax2(1) = subplot(1,2,1)
    for i = 1 : NSteps
        %Plot MotorAngle positive and negative for all steps together for each
        %cycle
        plot(T,MotorAngle_p{i,cycle} , 'Color', c(i,:) );hold on
        plot(T,MotorAngle_n{i,cycle} , 'Color', c(i,:) );
    end
    title([{Exp};{'MotorAngle'};param],'interpreter','none')
    
    %Plot MotorAngle positive and negative for all steps together for each
    %cycle
    
    ax2(2) = subplot(1,2,2)
    
    for i = 1 : NSteps
        % First, average the offset between -1 y -2 seconds
        offset_p{i,cycle} = nanmean( TailAngle_p{i,cycle}( find( and(T < -1,T > -2) ) ) );
        offset_n{i,cycle} = nanmean( TailAngle_n{i,cycle}( find( and(T < -1,T > -2) ) ) );
        
        %Subsctract the offset to put the the tail angle to zero
        TailAngle_p{i,cycle} = TailAngle_p{i,cycle} -offset_p{i,cycle};
        TailAngle_n{i,cycle} = TailAngle_n{i,cycle} -offset_n{i,cycle};
        
        plot(T,TailAngle_p{i,cycle} , 'Color', c(i,:));hold on
        plot(T,TailAngle_n{i,cycle} , 'Color', c(i,:) );
        ylim([-20 20 ])
    end
    title(['TailAngle cycle' int2str(cycle)])
    linkaxes(ax2,'x')
    
end
% So far, we have the figure for the motor angle and tail angle of 6 steps from 2 to 12 degrees
%% Plot in one figure all cycles (9) for every step separated (8 dif plots)

NSteps = 15;%change the  number according to protocol

Ncycles = 8%change the  number according to protocol


c = jet(NSteps); % color jet scale for the different steps

Step_fig = figure;
ax2(1) = subplot(4,4,1)
for i = 1 : NSteps
    plot(T,MotorAngle_p{i,Ncycles} , 'Color', c(i,:), 'LineWidth', 1);hold on
    plot(T,MotorAngle_n{i,Ncycles} , 'Color', c(i,:), 'LineWidth', 1 );
    
%     legend_info{2*i-1} = [int2str(2*i) ' deg']
    
    ylabel ('Motor Angle')
    xlim ([-30 30])
end
% l = legend(legend_info{:})
title([{Exp};{'MotorAngle'}; 'Navg = ' int2str(NavgTail); param] ,'interpreter','none')



for cycle = 1:Ncycles%change the initial number according to protocol
    
    for i = 1 : NSteps
        ax2(2) = subplot(4,4,1+i)% average the offset between -1 y -2 seconds
        offset_p{i,cycle} = nanmean( TailAngle_p{i,cycle}( find( and(T < -1,T > -2) ) ) );
        offset_n{i,cycle} = nanmean( TailAngle_n{i,cycle}( find( and(T < -1,T > -2) ) ) );
        
        TailAngle_p{i,cycle} = TailAngle_p{i,cycle} -offset_p{i,cycle};
        TailAngle_n{i,cycle} = TailAngle_n{i,cycle} -offset_n{i,cycle};
        
        plot(T,TailAngle_p{i,cycle} , 'Color', c(i,:), 'LineWidth', 1.5);hold on
        plot(T,TailAngle_n{i,cycle} , 'Color', c(i,:), 'LineWidth', 1.5 );
        ylim([-20 20 ])
        ylabel('Tail Angle')
        xlim ([-30 30])
        if i >= 1
        title([int2str(2*i) ' deg'])
        end

    end
    xlim ([-31 30])
end
xlim ([-30 30])
xlabel ('Time (sec)')
linkaxes(ax2,'x')
%% save workspace
% cd([path folder1])
% save('StepAnalysis.mat')
% savefig(findStep_fig,'findStep_fig.fig')
savefig(Step_fig,'StepAnalysis.fig')
% saveas(findStep_fig,'findStep','pdf')
saveas(Step_fig,'StepAnalysis','pdf')


%% Average along cycles

for i=1:NSteps
    for cycle=1:7
        tempMn(cycle,:)=MotorAngle_n{i,cycle};
        tempMp1(cycle,:)=MotorAngle_p{i,cycle};
        tempTn1(cycle,:)=TailAngle_n{i,cycle};
        tempTp1(cycle,:)=TailAngle_p{i,cycle};
    end
    MotorAngle_avg_n{i}=nanmean(tempMn);
    MotorAngle_avg_p{i}=nanmean(tempMp1);
    TailAngle_avg_n{i}=nanmean(tempTn1);
    TailAngle_avg_p{i}=nanmean(tempTp1);
    TailAngle_std_p{i}=nanstd(tempTp1);
    TailAngle_std_n{i}=nanstd(tempTn1);
    TailAngle_sem_p{i}=nanstd(tempTp1) / sqrt(size(TailAngle_std_p,2));
    TailAngle_sem_n{i}=nanstd(tempTn1) / sqrt(size(TailAngle_std_p,2));
    
end

% length(TailAngle_avg_p{8})
% t = 3000
% figure;

for i = 1:size(TailAngle_std_p,2)
errorbar(MotorAngle_avg_p{i}(t),TailAngle_avg_p{i}(t),TailAngle_std_p{i}(t));hold on
end
%% Plot average with SEM in bounded line

figure
c = jet (NSteps)
for i = 1:NSteps

boundedline(T,TailAngle_avg_p{i},TailAngle_sem_p{i}, 'cmap', c(i,:), 'alpha');
legend_info{(2*i-1)} = [num2str(2*i) ' deg'];
legend_info{(2*i)} = [num2str(2*i) ' deg'];

end
legend (legend_info)

%%
%% Display Averages

% figure
% for i = 1:NSteps
%     subplot(2,4,i)
%     plot(T, MotorAngle_avg_p{i})
%     hold on; plot(T, TailAngle_avg_p{i})
%     plot (T, TailAngle_avg_n{i})
%     plot(T, MotorAngle_avg_n{i})
%
%     title([int2str(2*i) ' deg'])
%     xlim([-30 30])
%     ylim([-20 0])
% end

%% Plot all cycles and their average for every step separated

AllCyclesAvg = figure
c = cool(Ncycles)
for i = 1:NSteps;
    
    
    subplot(5,3,i)
%     boundedline(T,TailAngle_avg_p{i},TailAngle_sem_p{i});hold on
%     boundedline(T,TailAngle_avg_n{i},TailAngle_sem_n{i});hold on

    for k = 2:Ncycles
        MTailAngle_p(k,:)= vec2mat([TailAngle_p{i,k}],4401);
        MTailAngle_n(k,:)= vec2mat([TailAngle_n{i,k}],4401);
        plot(T, MTailAngle_p(k,:)', 'Color', c(k,:));hold on
        plot(T, MTailAngle_n(k,:)', 'Color', c(k,:));
        legendInfo{(2*k-3)} = ['cycle n. ' num2str(k)];
        legendInfo{(2*k-2)} = ' ';
        
        
    end
    if i==1
        title ([{Exp};{'MotorAngle'}; 'Navg = ' int2str(NavgTail); param],'interpreter','none')
    end
    hold on
    plot(T,TailAngle_avg_p{:,i}, '-k', 'LineWidth',2);
    plot(T,TailAngle_avg_n{:,i}, '-k', 'LineWidth',2);
    

    if i > 1
        title([int2str(2*i) ' deg'])
    end
    ylim([-20 20])
    xlim([-30 30])
   
end

legend(legendInfo,'mean');
%%

cd([path folder1])
savefig(AllCyclesAvg,'AllCyclesAvg.fig')
% saveas(findStep_fig,'findStep','pdf')
saveas(AllCyclesAvg,'AllCyclesAvg','pdf')

%%

clear all
close all



%%


MTailAngle_avg_p = vec2mat([TailAngle_avg_p{:}],4401);
MTailAngle_std_p = vec2mat([TailAngle_std_p{:}],4401);
MTailAngle_avg_n = vec2mat([TailAngle_avg_n{:}],4401);
MTailAngle_std_n = vec2mat([TailAngle_std_n{:}],4401);

MMotorAngle_avg_p = vec2mat([MotorAngle_avg_p{:}],4401);

MMotorAngle_avg_n = vec2mat([MotorAngle_avg_n{:}],4401);

ErrorBar = figure
i = 2500;
errorbar(MMotorAngle_avg_p(:,i),MTailAngle_avg_p(:,i),MTailAngle_std_p(:,i), 's','MarkerSize', 8, 'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(MMotorAngle_avg_p(:,i),MTailAngle_avg_n(:,i),MTailAngle_std_n(:,i), 's','MarkerSize', 8, 'MarkerEdgeColor','cyan','MarkerFaceColor','cyan');

f = 4200
errorbar(MMotorAngle_avg_p(:,f),MTailAngle_avg_p(:,f),MTailAngle_std_p(:,f),'s','MarkerSize', 8, 'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(MMotorAngle_avg_p(:,f),MTailAngle_avg_n(:,f),MTailAngle_std_n(:,f), 's','MarkerSize', 8, 'MarkerEdgeColor','magenta','MarkerFaceColor','magenta');

xlabel ('Motor Angle (deg)')
ylabel ('Tail Angle (deg)')
legend ('left stim, early','right stim, early', 'left stim, late', 'right stim, late')

    title([{Exp};{'MotorAngle'}; param],'interpreter','none')

%%
cd([path folder1])
savefig(ErrorBar,'ErrorBar.fig')
% saveas(findStep_fig,'findStep','pdf')
saveas(ErrorBar,'ErrorBar','pdf')

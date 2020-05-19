%%
clear 

%%

close all
%% Open Folder and loading Data

path = '/home/ljp/Science/Projects/Rolling/2018-03/2018-03-20/';
folder = 'Run.01'
% 
% path = '/home/ljp/Science/Projects/Rolling/Data_balancoir/2018-01/'
% folder = 'ramp_90_500_4dpf_01'
load ([path folder '/data.mat'])

sep = find(path == filesep); %find the separations in the directory made by a /
Exp = [path(sep(end-2):end) folder]; % defines the name of the folder or directory to name the

cd (fullfile(path, folder))

v = VideoReader('video.avi');

Tf= v.NumberOfFrames;




%% Show one single frames without ellipses
figure %open a new figure

i = 180%fram index
%im = (imread(filelist(i).name));

im = read(v,i);

imshow(im)
%% Extract the angle of the 2 ellipses for the image sequence 

Tf= v.NumberOfFrames;%total number of frames
fi = 1;%initital frame
inc = 1;%increament in the number of frames, if you want to analyze all the sequence it should be = 1
ff = Tf;% final frame, it corresponds to the total lenght of the file list

tic
for i = fi : inc : ff
    if mod(i,100) == 0 %to display the number of frame analyzed in the command window
    disp(i);
    end
    im =read(v,i);

%     imgray = rgb2gray(im);
    [Bin]= ExtractTailH(im);% define [Bin] as the function "ExtractTailH" on image (im)
%     imshow(Bin) %show bin image on the figure
    hold on
    Es(i) = get_ellipses(Bin,2);
    Es_tot = get_ellipse(Bin);
%     draw_ellipse(Es.E1,[1 0 1]);%to draw ellipse E1 in RGB code color
%     pause (0.1)
%     draw_ellipse(Es.E2,[1 0 1]);;%to draw ellipse E2 in RGB code color
%     pause(0.1)
    Theta1(i) = Es(i).E1.theta*180/pi;  %call ellipse 1 angle in radians, E1 = head, and then transform to degrees
    Theta2 (i)= Es(i).E2.theta*180/pi; % call ellipse 2 angle in radians E2 = tail, and then transform to degrees
    angle_dif(i) = Theta1(i) - Theta2(i) ;%obtain delta between both angles in degrees


end
toc


%% Save the ellipses and angle data
path2 = '/home/ljp/Science/Projects/Rolling/Analysis';


mkdir ([path2 Exp]); %create a new directory to put the analysis in the folder ...Rolling/Analysis

cd([path2 Exp]);%change directory to the new created folder
save('angle_dif');%save data obtained from the subrstraction of the angle between the 2 ellipses





%% Plot angle delta and motor angle vs time stamp

Time = Data.TimeStamp; %rename all the variables contained in the 'Data' structure
TailH = Data.TailAngle - Data.TailAngle(1);% Tail H is the tail angle obtained directly in the program by Hugo's program
Tail2e = angle_dif - angle_dif(1); % Tail angle obtained using the 2 ellipses correction, minus the outset
N = 10 % number of data to average with the moving average method of smothing
Tail2e_smooth = smoothdata(Tail2e,'movmean',N);
Motor = Data.MotorAngle;


fig = figure
hold on;

yyaxis left
plot(Time (fi:inc:ff), TailH(fi:inc:ff), 'color', 'c')

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
legend('Tail H',['Tail 2e smoothen movmean' num2str(N)],'Motor Angle', 'Location','Southoutside')
title ([Exp])


%% Save figure as .fig and .pdf, save workspace
savefig ('Tail_vs_Motor.fig')
saveas (fig, 'Tail_vs_Motor.pdf')

save('workspace')



%% Show images and Draw ellipses on different frames. 

%% Draw the two ellipses in one single frame

figure %open a new figure
for i = 1 :10:100
    im = read (v,i);
    imgray = rgb2gray(im);
    [Bin]= ExtractTailH(imgray);% define [Bin] as the function "ExtractTailH" on image (im)
    imshow(Bin) %show bin image on the figure
    hold on
    Binrot = imrotate(Bin,20);
    Es = get_ellipses(Bin,2);
    draw_ellipse(Es.E1,[1 0 1]);
    draw_ellipse(Es.E2,[1 0 1]);
    pause(1)
end





%% Binarize, get and draw ellipse on one single frame


figure %open a new figure
i = 1
im = read(v,i);
imgray = rgb2gray(im)
[Bin]= ExtractTailH(imgray)% define [Bin] as the function "ExtractTailH" on image (im)
imshow(Bin) %show bin image on the figure
pause (0.5)
hold on
E = get_ellipse(Bin)
draw_ellipse(E,[1 0 1])




%%  Make an image sequence to show one ellipse
figure
for i = 1 : 10 : 10 %number of image file you want to open
    im = read(v,i);%convert image RGB to grayscale
    imshow(im) %show image on figure
    
    hold on %keep the image on the open window
    imgray = rgb2gray(im)
    [Bin]= ExtractTailH(imgray) %define [Bin] as the function writen on the script ExtractTailH.m over the im vector
    imshow(Bin)
    E (i)= get_ellipse(Bin);    % get ellipse from the treated image
    draw_ellipse(E(i),[1 0 0]); % draw elipse in red color RGB = [red= 1 green= 0 blue= 0]
    pause(1) %make a pause of 1 sec between images

end


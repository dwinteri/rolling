%%
clear all
close all
%% Open Folder and loading Data

path = '/home/ljp/Science/Projects/Natalia/Rolling/2017-12/2017-12-14/'

folder = 'run.1'

load ([path folder '/data.mat'])

cd (fullfile(path, folder))

v = VideoReader('video.avi');

Tf= v.NumberOfFrames;

%load data from the videostep2_16_6dpf_lor_1_long/images
% folder1 = '/step2_16_6dpf_1_long/'
% load([path folder 'data.mat'], 'Data')
%% Show single frames
figure %open a new figure

i = 10 %fram index
%im = (imread(filelist(i).name));

im = read(v,i);

imshow(im)
%% if images are in RGB they must be changed to 8bit or grayscale

for i = 1:Tf
   
    imcolor = read(v,i)
   
    imgray = rgb2gray(imcolor)
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
    [Bin]= ExtractTailH(im) %define [Bin] as the function writen on the script ExtractTailH.m over the im vector
    imshow(Bin)
    E (i)= get_ellipse(Bin);    % get ellipse from the treated image
    draw_ellipse(E(i),[1 0 0]); % draw elipse in red color RGB = [red= 1 green= 0 blue= 0]
    pause(1) %make a pause of 1 sec between images

end


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
    imgray = rgb2gray(im);
    [Bin]= ExtractTailH(imgray);% define [Bin] as the function "ExtractTailH" on image (im)
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
save([path folder '/AngleDif.mat'],'angle_dif')
save([path folder '/Es.mat'],'Es')



%% Plot angle delta and motor angle vs time stamp

figure
hold on;
plot(Data.TimeStamp (fi:inc:ff), angle_dif(fi:inc:ff))
% plot(angleE1(fi:inc:ff) - 90)
plot(Data.TimeStamp (fi:inc:ff), (Data.TailAngle(fi:inc:ff)-Data.TailAngle(1)))
plot (Data.TimeStamp(fi:inc:ff), Data.MotorAngle(fi:inc:ff))
legend('Tail Angle E2','Tail Angle H','Motor Angle')
% plot(-Data.TailAngle + 90)
% savefig ('/media/ljp/800a0dce-e801-4fcc-83b9-3f25412261e5/Science/Projects/Natalia/rolling/january/4dpf/stair_90_4dpf_04/AngleDif.fig')





%% Draw the two ellipses in one single 
figure %open a new figure
for i = 1 :10:100
    im = imread(filelist(i).name);
    [Bin]= ExtractTailH(im);% define [Bin] as the function "ExtractTailH" on image (im)
    imshow(Bin) %show bin image on the figure
    hold on
    Binrot = imrotate(Bin,20);
    Es = get_ellipses(Bin,2);
    draw_ellipse(Es.E1,[1 0 1]);
    draw_ellipse(Es.E2,[1 0 1]);
    pause(1)
end

Es.E1.theta*180/pi %call ellipse 1 angle in radians, E1 = head


Es.E2.theta*180/pi % call ellipse 2 angle in radians E2 = tail




%%




% for 
% while
% 
% if
% 
% switch name
%     case name == 'Natalia'
% end
% 
% save 
% saveas
% 
% repmat
% 
% x = get(gco,'xdata')
% y = get(gco,'ydata')
% figure;plot(x,y)
% 
% title
% xlabel
% ylabel













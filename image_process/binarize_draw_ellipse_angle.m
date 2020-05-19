%%
clear all
close all
%% Open Folder and loading Data

path = '//home/ljp/Science/Projects/Natalia/Rolling/2018-03/2018-03-20/'
folder = 'Run.03'
cd ('/home/ljp/Science/Projects/Natalia/Rolling/2018-03/2018-03-20/Run.03/images')
filelist = dir('*.tif')%show list of files on the current directory

%load data from the videostep2_16_6dpf_lor_1_long/images
% folder1 = '/step2_16_6dpf_1_long/'
% load([path folder 'data.mat'], 'Data')
%% Show single frames
figure %open a new figure
i = 10 %fram index
%im = (imread(filelist(i).name));
im = (imread(filelist(i).name));
imshow(im)
%% Binarize, get and draw ellipse on one single frame


figure %open a new figure
i = 780
im = (imread(filelist(i).name));
[Bin]= ExtractTailH(im)% define [Bin] as the function "ExtractTailH" on image (im)
imshow(Bin) %show bin image on the figure
pause (0.1)
hold on
E = get_ellipse(Bin)
draw_ellipse(E,[1 0 1])




%%  Make an image sequence to show one ellipse
figure
for i = 1500 : 300 : 4500 %number of image file you want to open
    im = ( imread(filelist(i).name) );%convert image RGB to grayscale
    imshow(im) %show image on figure
    hold on %keep the image on the open window
    [Bin]= ExtractTailH(im) %define [Bin] as the function writen on the script ExtractTailH.m over the im vector
    imshow(Bin)
    E (i)= get_ellipse(Bin);    % get ellipse from the treated image
    draw_ellipse(E(i),[1 0 0]); % draw elipse in red color RGB = [red= 1 green= 0 blue= 0]
    pause(1) %make a pause of 1 sec between images

end


%% Extract the angle of the 2 ellipses for the image sequence 
Tf= length(filelist);%total number of frames
fi = 1;%initital frame
inc = 1;%increament in the number of frames, if you want to analyze all the sequence it should be = 1
ff = Tf;% final frame, it corresponds to the total lenght of the file list

tic
for i = fi : inc : ff
    if mod(i,100) == 0 %to display the number of frame analyzed in the command window
    disp(i);
    end
    im = imread(filelist(i).name);
    [Bin]= ExtractTailH(im);% define [Bin] as the function "ExtractTailH" on image (im)
%     imshow(Bin) %show bin image on the figure
    hold on
    Es(i) = get_ellipses(Bin,2);
    Es_tot = get_ellipse(Bin);
%     draw_ellipse(Es.E1,[1 0 1]);
%     pause (0.1)
%     draw_ellipse(Es.E2,[1 0 1]);
%     pause(0.1)
    Theta1(i) = Es(i).E1.theta*180/pi;  %call ellipse 1 angle in radians, E1 = head
    Theta2 (i)= Es(i).E2.theta*180/pi; % call ellipse 2 angle in radians E2 = tail
    angle_dif(i) = Theta1(i) - Theta2(i) ;%obtain delta between both angles in degrees
%     angleE1(i) = Es_tot.theta*180/pi  ; %obtain angle of one ellipse (get.ellipse)

end
toc


%% Save the ellipses and angle data
save([path folder 'AngleDif.mat'],'angle_dif')
save([path folder 'Es.mat'],'Es')



%% Plot angle delta and motor angle vs time stamp

figure
hold on;
plot(Data.TimeStamp (fi:inc:ff), angle_dif(fi:inc:ff)-4)
% plot(angleE1(fi:inc:ff) - 90)

plot (Data.TimeStamp(fi:inc:ff), Data.MotorAngle(fi:inc:ff))
legend('Tail Angle E2','Motor Angle')
% plot(-Data.TailAngle + 90)
% savefig ('/media/ljp/800a0dce-e801-4fcc-83b9-3f25412261e5/Science/Projects/Natalia/rolling/january/4dpf/stair_90_4dpf_04/AngleDif.fig')





%% Obtain angle difference of two ellipses

%call Es to show ellipses properties



Es.E1.theta*180/pi %call ellipse 1 angle in radians, E1 = head


Es.E2.theta*180/pi % call ellipse 2 angle in radians E2 = tail


Es.E2.theta*180/pi - Es.E1.theta*180/pi 

%% Draw the two ellipses in one single 
figure %open a new figure
for i = 1200 :300:4800
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


%     0.0856    0.1778    0.4112 





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













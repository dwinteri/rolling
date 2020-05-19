%% Binarize and get and draw ellipse on an image sequence

%open file
cd('/home/ljp/Natalia/rolling/december/ramp_stim/ramp/ramp_500_8dpf_06/Images')%cd = change directory
filelist = dir('*.tif');%show list of files on the current directory
figure %open a new figure


%open a image sequence
figure
for i = 1 : 10 : 100 %number of image file you want to open
    im = rgb2gray( imread(filelist(i).name) );%convert image RGB to grayscale
    imshow(im) %show image on figure
    hold on %keep the image on the open window
    [Treat]=ExtractTailH(im) %define [Treat] as the function writen on the script ExtractTailH.m over the im vector
    E = get_ellipse(Treat);    % get ellipse from the treated image
    draw_ellipse(E,[1 0 0]); % draw elipse in red color RGB = [red= 1 green= 0 blue= 0]
    pause(1) %make a pause of 1 sec between images

end

%% Binarize and  get and draw ellipse on one image

cd('/home/ljp/Natalia/rolling/december/ramp_stim/ramp/ramp_500_8dpf_06/Images')%cd = change directory
filelist = dir('*.tif');%show list of files on the current directory

figure %open a new figure

%open a image sequence
    i = 1  %number of image file you want to open
    im = rgb2gray( imread(filelist(i).name) );%convert image RGB to grayscale
    imshow(im) %show image on figure
    hold on %keep the image on the open window
    [Treat]=ExtractTailH(im) %
    E = get_ellipse(Treat);    % get ellipse from the treated image
    draw_ellipse(E,[0 0 1]); % draw elipse in red color RGB = [red= 1 green= 0 blue= 0]
    pause(1) %make  pause of 1 sec between images
%% Binarize and get and draw ellipse on one binarized image

cd('/home/ljp/Natalia/rolling/January_2018/ramp_80_500_8dpf_04/images')%cd = change directory
filelist = dir('*.tif')

%%
figure %open a new figure
i = 1
im = rgb2gray(imread(filelist(i).name));
[Bin]= ExtractTailH(im)% define [Bin] as the function "ExtractTailH" on image (im)
imshow(Bin) %show bin image on the figure
hold on
E = get_ellipse(Bin)
draw_ellipse(E,[1 0 1])
pause(1)

%% Binarize and  get and draw 2 ellipses on one image

cd('/home/ljp/Natalia/rolling/december/ramp_stim/ramp/ramp_500_8dpf_06/Images')%cd = change directory
filelist = dir('*.tif')
%%
clear angleE2 angleE1
figure %open a new figure
inc = 50;
for i = 1 :inc:9000
    clear Es Es_tot 
    im = rgb2gray(imread(filelist(i).name));
    [Bin]= ExtractTailH(im);% define [Bin] as the function "ExtractTailH" on image (im)
%    imshow(Bin) %show bin image on the figure
    hold on
    Es = get_ellipses(Bin,2);
    Es_tot = get_ellipse(Bin);
%     draw_ellipse(Es.E1,[1 0 1]);
%     draw_ellipse(Es.E2,[1 0 1]);
%     pause(1)
    angleE2(i) = Es.E2.theta*180/pi - Es.E1.theta*180/pi ;
    angleE1(i) = Es_tot.theta*180/pi  ;

end
figure;plot(angleE2(1:inc:end));hold on;plot(angleE1(1:inc:end))
legend('E2','E1')

%% Obtain angle difference of two ellipses

%call Es to show ellipses properties



Es.E1.theta*180/pi %call ellipse 1 angle in radians, E1 = head


Es.E2.theta*180/pi % call ellipse 2 angle in radians E2 = tail


Es.E2.theta*180/pi - Es.E1.theta*180/pi 

%%
figure %open a new figure
for i = 1 
    im = rgb2gray(imread(filelist(i).name));
    [Bin]= ExtractTailH(im);% define [Bin] as the function "ExtractTailH" on image (im)
    imshow(Binrot) %show bin image on the figure
    hold on
    Binrot = imrotate(Bin,20);
    Es = get_ellipses(Binrot,2);
    draw_ellipse(Es.E1,[1 0 1]);
    draw_ellipse(Es.E2,[1 0 1]);
    pause(1)
end

Es.E1.theta*180/pi %call ellipse 1 angle in radians, E1 = head


Es.E2.theta*180/pi % call ellipse 2 angle in radians E2 = tail
























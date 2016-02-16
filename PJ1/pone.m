clc;
clear all;
close all;

 
srcFiles = dir([pwd '\EnterExitCrossingPaths2cor\*.jpg']);
        % the folder in which ur images exists
   
thresh = 15;
    
for i = 60 %1 : 400
    
    filename1 = strcat(pwd, '\EnterExitCrossingPaths2cor\', srcFiles(i).name);
    I1 = imread(filename1);         
    bg_bw=rgb2gray(I1);
    
    fr_size = size(bg_bw);
    width = fr_size(2);
    height = fr_size(1);
    fg = zeros(height, width);
    
    filename2 = strcat(pwd, '\EnterExitCrossingPaths2cor\', srcFiles(i+1).name);
    I2 = imread(filename2);
    f2_bw = rgb2gray(I2);
    
    filename3 = strcat(pwd, '\EnterExitCrossingPaths2cor\', srcFiles(i+2).name);
    I3 = imread(filename3);
    %figure, imshow(I);
              
    fr_bw = rgb2gray(I3); 
    %figure, imshow(Igrey);
    
    fr_diff = (abs(double(fr_bw) - double(bg_bw)))/2;  % Cast operands as double to avoid negative overflow
       
    Mask = image_threshold ( fr_diff, thresh );
    
   figure(1), 
   % subplot(1,2,1), imshow(I2)
   % subplot(2,2,2),imshow(fr_bw)
   % subplot(2,2,3),imshow(uint8(fr_diff))
   % subplot(1,2,2),
    imshow(Mask) 
end
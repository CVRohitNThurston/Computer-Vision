%% 1D Gaussian filter

clc;
clear all;
close all;

 
srcFiles = dir([pwd '\EnterExitCrossingPaths2cor\*.jpg']);
% the folder in which ur images exists

%% prep display
figure(1), placeplot(1,1); fax(1) = gca;
figure(2), placeplot(2,2); fax(2) = gca;
figure(3), placeplot(3,3); fax(3) = gca;
figure(4), placeplot(4,4); fax(4) = gca;

%% gaussian kernel
sigma = 1;

% Determine filter length
filterLength = ceil(5*(sigma)) + mod(ceil(5*(sigma))-1,2);
n = (filterLength - 1)/2;
x = -n:n;

% Create 1-D Gaussian Kernel
c = 1/(sqrt(2*pi)*sigma);
gaussKernel = c * exp(-(x.^2)/(2*sigma^2));

% Normalize to ensure kernel sums to one
gaussKernel = gaussKernel/sum(gaussKernel);

% Create 1-D Derivative of Gaussian Kernel
derivGaussKernel = gradient(gaussKernel);
derivGaussKernel = derivGaussKernel/sum(abs(derivGaussKernel));
        
thresh = 10;

numFrames = length(derivGaussKernel);
%% process    
for i = 1 : 400
    
    for f = 1:numFrames
        filename = strcat([pwd '\EnterExitCrossingPaths2cor\'],srcFiles(i+f-1).name); 
        I = rgb2gray(imread(filename));
        bg(:,:,f) = double(I); % image-time matrix: (row,column,frame)
    end
       
    fr_size = size(bg);
    width = fr_size(2);
    height = fr_size(1);
    fg = zeros(height, width);
    
    %% 2D smoothing
    smoother = 1;
    if (smoother)
       % choose filter 
       smooth_type = 2;
       switch(smooth_type)
           case 1, % 3x3 box filter
               filt = ones([3 3]);
           case 2, % 5x5 box filter
               filt = ones([5 5]);
           case 3, % 2D Gaussian
               filt = bsxfun(@times,gaussKernel,gaussKernel.');
       end
       %normalize
       filt = filt ./ sum(sum(abs(filt)));
       
       % frame by frame convolution
       bgprime = zeros(size(bg));
       for f = 1:numFrames
           bgprime(:,:,f) = conv2(bg(:,:,f),filt,'same');
       end
    else
        bgprime = bg;
    end
    
    
   %% Convolve with 1D gaussian in the temporal domain
   frameFactor = bsxfun(@times,double(bgprime),shiftdim(derivGaussKernel,-1));
   fr_diff = sum(frameFactor,3);
   
   %% Threshold
   Mask = image_threshold ( fr_diff, thresh );
     
   % subplot(1,2,1), imshow(I2)
   % subplot(2,2,2), imshow(fr_bw)
   % subplot(2,2,3),
   % subplot(1,2,2),
   
   %% Display
   imshow(uint8(bg(:,:,ceil(numFrames/2))),'Parent',fax(1));
   imshow(uint8(bgprime(:,:,ceil(numFrames/2))),'Parent',fax(2));
   imshow(uint8(fr_diff),'Parent',fax(3));
   imshow(Mask,'Parent',fax(4))
   pause(0.05);
end
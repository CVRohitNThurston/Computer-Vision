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
tsigma = 4.2;

% Determine filter length
filterLength = ceil(5*(tsigma)) + mod(ceil(5*(tsigma))-1,2);
n = (filterLength - 1)/2;
x = -n:n;

% Create 1-D Gaussian Kernel
c = 1/(sqrt(2*pi)*tsigma);
gaussKernel = c * exp(-(x.^2)/(2*tsigma^2));

% Normalize to ensure kernel sums to one
gaussKernel = gaussKernel/sum(gaussKernel);

%% Create 1-D Derivative of Gaussian Kernel
derivGaussKernel = gradient(gaussKernel);
derivGaussKernel = derivGaussKernel/sum(abs(derivGaussKernel));
        

%% choose time domain kernel
kerneltype = 2;
switch(kerneltype)
    case 1,
        timekernel = 0.5 .* [-1 0 1];
    case 2,
        timekernel = derivGaussKernel;
end

numFrames = length(timekernel);

%% declare spatial Gaussain kernel
ssigma = 4.2;

% Determine filter length
filterLength = ceil(5*(ssigma)) + mod(ceil(5*(ssigma))-1,2);
n = (filterLength - 1)/2;
x = -n:n;

% Create 1-D Gaussian Kernel
c = 1/(sqrt(2*pi)*ssigma);
gaussKernel = c * exp(-(x.^2)/(2*ssigma^2));

% Normalize to ensure kernel sums to one
spaceKernel = gaussKernel/sum(gaussKernel);


%% process    
for i = 60 % 1 : 400
    
    for f = (0-floor(numFrames/2)):(0+floor(numFrames/2))
        filename = strcat([pwd '\EnterExitCrossingPaths2cor\'],srcFiles(i+f-1).name); 
        I = rgb2gray(imread(filename));
        bg(:,:,f+1+floor(numFrames/2)) = double(I); % image-time matrix: (row,column,frame)
    end
       
    fg = zeros(size(bg));
    
    %% 2D smoothing filter
    % choose filter
    smooth_type = 1;
    switch(smooth_type)
        case 1, % no smoothing
            filt = 1;
        case 2, % 3x3 box filter
            filt = ones([3 3]);
        case 3, % 5x5 box filter
            filt = ones([5 5]);
        case 4, % 2D Gaussian
            filt = bsxfun(@times,spaceKernel,spaceKernel.');
    end
    
    %normalize
    filt = filt ./ sum(sum(abs(filt)));
    
    % frame by frame convolution
    bg_smooth = zeros(size(bg));
    for f = 1:numFrames
        bg_smooth(:,:,f) = conv2(bg(:,:,f),filt,'same');
    end    
    
    
   %% Correlate with 1D gaussian in the temporal domain
   frameFactor = bsxfun(@times,double(bg_smooth),shiftdim(timekernel,-1));
   fr_diff = abs(sum(frameFactor,3));

   
   %% Get noise estimate for thresholding
   % get variance per pixel
   bg_var = sum( (bsxfun(@minus,bg_smooth,mean(bg_smooth,3)).^2), 3)...
       ./ (numFrames - 1);
   % choose median variance across the camera : assumes that less than half
   % of the pixels are confounded with motion
   bg_thresh = 5 * sqrt( median( reshape(bg_var,1,[]) ));
   
   bg_thresh = 6;
   
   %% Threshold
   Mask = image_threshold ( fr_diff, bg_thresh );

   
   %% Display
   imshow(uint8(bg(:,:,ceil(numFrames/2))),'Parent',fax(1));
   imshow(uint8(bg_smooth(:,:,ceil(numFrames/2))),'Parent',fax(2));
   imshow(uint8(fr_diff)*30,'Parent',fax(3));
   imshow(Mask,'Parent',fax(4))
   pause(0.0000001);
end
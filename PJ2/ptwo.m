clc;
clear all;
close all;

srcFiles = dir([pwd '\DanaOffice\*.jpg']);
% srcFiles = dir([pwd '\DanaHallWay2\*.jpg']);
numFrames = size(srcFiles);
ImGreySet= zeros(340,512,10,'uint8');
for f = 1:numFrames(1)
    filename = strcat([pwd '\DanaOffice\'],srcFiles(f).name);
%     filename = strcat([pwd '\DanaHallWay2\'],srcFiles(f).name);
    I(:,:,:,f) = imread(filename);
    ImGreySet(:,:,f) = uint8(rgb2gray(I(:,:,:,f))); % image-time matrix: (row,column,frame)
end

index = 3;
I1 = I(:,:,:,1+index);
I2 = I(:,:,:,2+index);

Igrey1 = ImGreySet(:,:,1+index);
Igrey2 = ImGreySet(:,:,2+index);

%% Apply Corner Detector for the image Set
warning off;
CornerSet1 = harris(Igrey1,1,4,25000,0);
CornerSet2 = harris(Igrey2,1,4,25000,0);
warning on;

[r, c]=size(Igrey1);
Thresh = 0.8;

k = int8(1);

%% NCC
for i = 1:length(CornerSet1)
    % Choosing a neighborhood for the corner point of first image
    CornerPoint1 = CornerSet1(i,:);
    if (CornerPoint1(1)<=10) || (CornerPoint1(2)<= 10) || (CornerPoint1(1)> r-10) || (CornerPoint1(2)> c-10)
        continue;
    end
    nbhd1 = I1((CornerPoint1(1)-10):(CornerPoint1(1)+10),(CornerPoint1(2)-10):(CornerPoint1(2)+10));
    
    NCCArray = zeros(1,length(CornerSet2));
    for j = 1: length(CornerSet2)
        % Choosing a neighborhood for the corner point of first image
        CornerPoint2 = CornerSet2(j,:);
        if ((CornerPoint2(1)<11) || (CornerPoint2(2)< 11) || (CornerPoint2(1)> r-10) || (CornerPoint2(2)> c-10) )
            continue;
        end
        nbhd2 = I2((CornerPoint2(1)-10):(CornerPoint2(1)+10),(CornerPoint2(2)-10):(CornerPoint2(2)+10));
        NCC = normxcorr2(nbhd1, nbhd2);
        NCCArray(1,j) = NCC(21,21);
    end
    [LargestNCC, jIndex]= max(NCCArray(:));
    if  LargestNCC > Thresh
        CorrespMap(k,1:2) = [i jIndex];
        k=k+1;
    end
end
Np = length(CorrespMap);
%% Plot initial point correspondences
Cset1Index = CorrespMap(:,1);  Cset2Index = CorrespMap(:,2);

figure;
CombinedImage = cat(2,I1,I2);
imshow(CombinedImage);
hold on

xs1 = CornerSet1(:,2); 
ys1 = CornerSet1(:,1); 
xs2 = CornerSet2(:,2); 
ys2 = CornerSet2(:,1);

plot(xs1, ys1, 'gs','LineWidth',2);
plot(xs2+c+1, ys2, 'rx','LineWidth',2);
title('Initial Point Correspondences')

%% Plotting initial point correspondences
hold on;
for i = 1: Np
    
%     if(last_inliers(i))
    %
    x1(i) = CornerSet1(Cset1Index(i),2);
    y1(i) = CornerSet1(Cset1Index(i),1);
    
    %
    x2(i) = CornerSet2(Cset2Index(i),2);
    y2(i) = CornerSet2(Cset2Index(i),1);

    plot([x1(i) x2(i)+c+1],[y1(i) y2(i)]);
%     end
end

%% RANSAC to find Homography

h_ransac = ones([8 1]);
last_error = inf;
last_inliers = 0;
while(sum(last_inliers) < 8)
    for jj = 1:Np
        
        % choose 4 random points at a time
        i = randi([1 Np], [1 4]);
        
        x = x1(i);
        y = y1(i);
        xp = x2(i);
        yp = y2(i);
        
        % build A matrix
        for ii = 4:-1:1
            A(2*ii-1:2*ii,:) = [x(ii) y(ii) 1 0 0 0 -x(ii)*xp(ii) -y(ii)*xp(ii);...
                0 0    0 x(ii) y(ii) 1 -x(ii)*yp(ii) -y(ii)*yp(ii)];
        end
        
        % build b matrix
        b = reshape([xp;yp],8,[]);
        
        % get homography estimate
        warning off;
        h_est = A \ b;
        warning on;
        h_est = reshape([h_est; 1],3,3).';
        
        %     % show final image
        %     stitched_im = stitch_images(Igrey1,Igrey2,h_est);
        
        % get error of the estimate
        thresh = 5;
        [h_err, inliers] = h_error(h_est,x1,y2,x2,y2,thresh);
        
        % update estimate
        if(h_err < last_error)
            last_h_est = h_est;
            last_error = h_err;
            last_inliers = inliers;
        end
    end
end

%% get final estimate

% choose inliers
x = x1(last_inliers);
y = y1(last_inliers);
xp = x2(last_inliers);
yp = y2(last_inliers);

% build A matrix
for ii = length(xp):-1:1
    A(2*ii-1:2*ii,:) = [x(ii) y(ii) 1 0 0 0 -x(ii)*xp(ii) -y(ii)*xp(ii);...
        0 0    0 x(ii) y(ii) 1 -x(ii)*yp(ii) -y(ii)*yp(ii)];
end

% build b matrix
b = reshape([xp;yp],length(xp)*2,[]);

% get final homography estimate
h_est = A \ b;
h_ransac = reshape([h_est; 1],3,3).';

last_h_ransac = last_h_est;
%% MSAC to check

[tform, inlierPoints1, inlierPoints2] = estimateGeometricTransform([x1(:), y1(:)],[x2(:), y2(:)],'projective');
h_msac = (tform.T).';

%% Inliers(green) found by MSAC
figure;
imshowpair(Igrey1(:,:), Igrey2(:,:),'montage');
hold on;
plot(x1(:), y1(:), 'rx', 'LineWidth', 2);
plot(x2(:)+c+1, y2(:), 'bx', 'LineWidth', 2);
scatter(inlierPoints1(:,1), inlierPoints1(:, 2), 18,'g', 'fill','LineWidth', 5);
scatter(inlierPoints2(:,1)+c+1, inlierPoints2(:, 2),18, 'g','fill', 'LineWidth', 5);
title('Inlier Points Determined by M-SAC');


%% Inliers(green) found by RANSAC
figure;
imshowpair(Igrey1(:,:), Igrey2(:,:),'montage');
hold on;
plot(x1(:), y1(:), 'rx', 'LineWidth', 2);
plot(x2(:)+c+1, y2(:), 'bx', 'LineWidth', 2);
hold on
for i = 1: sum(last_inliers)
    plot([x(i) xp(i)+c+1],[y(i) yp(i)]);
end
title('Inlier Points Determined by RANSAC');

%% Warp one image onto the other one
%% Stitch 
[stitched_im, IM] = stitch_images(Igrey1,Igrey2,h_ransac);
figure; imshow(uint8(IM.source)); title('Source Image');
figure; imshow(uint8(IM.dest)); title('Destination ');
figure; imshow(uint8(IM.stitched)); title('The final stitched image');


%% Initialize
% clc;
% clear all;
% close all;

numFrames = 2;
imdir = '\Cones';
% imdir = '\Cast';

I1 = imread([pwd imdir '\left.jpg']);
I2 = imread([pwd imdir '\right.jpg']);
Igrey1 = uint8(rgb2gray(I1));
Igrey2 = uint8(rgb2gray(I2));

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
        clearvars A;
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

%% get final Homography estimate 

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

%% RANSAC to find Fundamental Matrix

last_error = inf;
last_inliers = 0;
while(sum(last_inliers) < 8)
    for jj = 1:Np
        
        %% choose 8 random points at a time
        i = randi([1 Np], [1 8]);
            
        x = x1(i);
        y = y1(i);
        xp = x2(i);
        yp = y2(i);
        
        %% build normalization matrices
        mux = mean(x); muxp = mean(xp);
        muy = mean(y); muyp = mean(yp);
        sig = sum(sqrt((x-mux).^2+(y-muy).^2))/(length(x)*sqrt(2));
        sigp = sum(sqrt((xp-muxp).^2+(yp-muyp).^2))/(length(xp)*sqrt(2));
        
        T = [1/sig, 0, -mux/sig; 0 1/sig, -muy/sig; 0 0 1];
        Tp = [1/sigp, 0, -muxp/sigp; 0 1/sigp, -muyp/sigp; 0 0 1];
        
        % get normalized points
        p = permute([x;y;ones(size(x))],[1 3 2]);
        pp = permute([xp;yp;ones(size(xp))],[1 3 2]);
        
        for i = size(p,3):-1:1; phd(:,:,i) = T*p(:,:,i); end;
        for i = size(pp,3):-1:1; pphd(:,:,i) = Tp*pp(:,:,i); end;
        
        ph = permute(phd, [1 3 2]);
        pph = permute(pphd, [1 3 2]);
              
        %% build data matrix, A
        A = [ph(1,:).'.*pph(1,:).',...
            ph(1,:).'.*pph(2,:).',...
            ph(1,:).',...
            ph(2,:).'.*pph(1,:).',...
            ph(2,:).'.*pph(2,:).',...
            ph(2,:).',...
            pph(1,:).',...
            pph(2,:).',...
            ones([size(ph,2) 1])];

        %% get the entries of F-normalized from the final eigenvector of 
        [U, D, V] = svd(A);
        
        Fn = reshape(V(:,end),3,3).';
        
        % verify that Fn relates points well
        for i = size(ph,2)
           S(i) = [ph(1,i) ph(2,i) 1] * Fn * [pph(1,i); pph(2,i); 1]; 
        end
        disp(sum(S(i)));
        
        % force F-normalized to be singular
        [Uf, ~, Vf] = svd(Fn); 
        Df = svd(Fn); 
        Fp = Uf*diag([Df(1:end-1); 0])*Vf.';
        
        % verify that Fp relates points well
        for i = size(ph,2)
           S(i) = [ph(1,i) ph(2,i) 1] * Fp * [pph(1,i); pph(2,i); 1]; 
        end
        disp(sum(S(i)));
        
        % denormalize F
        F = (inv(Tp).') * Fp * (T);
          
        %% find epic-polar-lines
        for i = 1:8
            l2(:,i) = F * [x(i); y(i); 1];
            l1(:,i) = F.' * [xp(i); yp(i); 1];
        end
        
        %% 
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

%% Stitch 
[stitched_im, IM] = stitch_images(Igrey1,Igrey2,h_ransac);
figure; imshow(uint8(IM.source)); title('Source Image');
figure; imshow(uint8(IM.dest)); title('Destination ');
figure; imshow(uint8(IM.stitched)); title('The final stitched image');

%% Create disparity map




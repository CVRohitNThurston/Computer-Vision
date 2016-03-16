clc;
clear all;
close all;

% I2 = imread('NYC.jpg');
I2 = imread('Barcelona.jpg');
I1 = imread('trump.jpg');
figure;
imshow(I1);
figure;
imshow(I2);

Igrey1 = uint8(rgb2gray(I1));
Igrey2 = uint8(rgb2gray(I2));

%% Manual corner points input

[x1,y1] = input_cornerpts(Igrey1);
[x2,y2] = input_cornerpts(Igrey2);

x1 = x1(1:4);
x2 = x2(1:4);
y1 = y1(1:4);
y2 = y2(1:4);

%% Homography
A = zeros(size(x1,1)*2, 9);
for i = 1:size(x1,1)
A(2*i-1,:) = [x1(i),y1(i),1, 0,0,0, -x1(i)*x2(i),-x2(i)*y1(i),-x2(i)];
A(2*i,:) = [0,0,0, x1(i),y1(i),1, -x1(i)*y2(i),-y2(i)*y1(i),-y2(i)];
end
% 
% for ii = length(x2):-1:1
%     A(2*ii-1:2*ii,:) = [x1(ii) y1(ii) 1 0 0 0 -x1(ii)*x2(ii) -y1(ii)*x2(ii);...
%         0 0    0 x1(ii) y1(ii) 1 -x1(ii)*y2(ii) -y1(ii)*y2(ii)];
% end
% 
% % build b matrix
% b = reshape([x2;y2],length(x2)*2,[]);
% 
% H = A \ b;
% H = reshape([H; 1],3,3).';
[U,S,V] = svd(A);
h = V(:,9);
H = [h(1),h(2),h(3);h(4),h(5),h(6);h(7),h(8),h(9)];


%% Warp the images
stitched_im = stitch_images(Igrey1,Igrey2,H,100);
figure;
imshow(stitched_im);
title('The final stitched image');



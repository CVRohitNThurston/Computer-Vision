function [ stitched_Im, pix ] = stitch_images( Source_Im, Dest_Im, H)

%% find min/max x/y coordinates

% find boundaries
p1s = [1; 1; 1];
p2s = [size(Source_Im,2); 1; 1];
p3s = [1; size(Source_Im,1); 1];
p4s = [size(Source_Im,2); size(Source_Im,1); 1];

p1d = H*p1s; p1d = p1d/p1d(3);
p2d = H*p2s; p2d = p2d/p2d(3);
p3d = H*p3s; p3d = p3d/p3d(3);
p4d = H*p4s; p4d = p4d/p4d(3);

% get min/max coordinates in destination frame
minx = floor(min([p1s(1), p2s(1), p3s(1), p4s(1), p1d(1), p2d(1), p3d(1), p4d(1)]));
miny = floor(min([p1s(2), p2s(2), p3s(2), p4s(2), p1d(2), p2d(2), p3d(2), p4d(2)]));
maxx = ceil(max([p1s(1), p2s(1), p3s(1), p4s(1), p1d(1), p2d(1), p3d(1), p4d(1)]));
maxy = ceil(max([p1s(2), p2s(2), p3s(2), p4s(2), p1d(2), p2d(2), p3d(2), p4d(2)]));

% get coords for the destination frame
[xi, yi] = meshgrid(minx:maxx,miny:maxy);
h = inv(H);

% get corresponding coordinates for the source frame using the homography
xx = (h(1,1)*xi+h(1,2)*yi+h(1,3))./(h(3,1)*xi+h(3,2)*yi+h(3,3));
yy = (h(2,1)*xi+h(2,2)*yi+h(2,3))./(h(3,1)*xi+h(3,2)*yi+h(3,3));

% get boolean 2D array of which pixels are valid
% source
pixs = (zeros(size(xx)));
pixs = (xx > 1) & (xx < size(Source_Im,2)) & (yy > 1) & (yy < size(Source_Im,1));

    
%destination
pixd = zeros(size(xx));
xd = 1:size(Dest_Im,2);
yd = 1:size(Dest_Im,1);
ydd = yd - miny + 1;
xdd = xd - minx + 1;
pixd(ydd, xdd) = 1;


% feathering
sigma = 3;
feather = bsxfun(@times,...
    normpdf(ydd,mean(ydd),size(Dest_Im,1)/(2*sigma)).',...
    normpdf(xdd,mean(xdd),size(Dest_Im,2)/(2*sigma)));
pixd(ydd,xdd) = feather;

pixw = bsxfun(@times,...
    normpdf(xx,mean([1 size(Source_Im,2)]),size(Source_Im,2)/(2*sigma)),...
    normpdf(yy,mean([1 size(Source_Im,1)]),size(Source_Im,1)/(2*sigma)));
pixs = pixs .* pixw;

% interpolate source image values to map them to the destination image
Source_mapped = uint8(interp2(double(Source_Im),xx,yy));

% map the destination to the expanded coordinate system
Dest_mapped = uint8(zeros(size(Source_mapped)));
Dest_mapped(ydd, xdd) = Dest_Im;


% average the overlapping pixels
Source_mapped = double(Source_mapped);
Dest_mapped = double(Dest_mapped);
pixd = double(pixd);
pixs = double(pixs);

% for i=size(xi,1);
%     for j=size(xi,2);
%         stitched_Im(i,j) = ( ((Source_mapped(i,j)) + (Dest_mapped(i,j) )) / ((pixs(i,j)) + (pixd(i,j))) );
%     end
% end

sum_im = Source_mapped.*pixs + Dest_mapped.*pixd;
pix = pixs + pixd;
overlap_pix = pix > 0;
sum_im(overlap_pix) = sum_im(overlap_pix) ./ pix(overlap_pix);

stitched_Im = uint8(sum_im);

% display
figure; imshow(uint8(Source_mapped))
figure; imshow(uint8(Dest_mapped));
% figure; imshow(stitched_Im);


end


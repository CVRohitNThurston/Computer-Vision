function [ psr ] = PSR( filterResponse )
% calculate the PSR (Peak to Sidelobe Ratio)
maxresponse = max(filterResponse(:));
windowSize = 11;
 % maximum value of filter response
[xx, yy] = find(filterResponse == maxresponse, 1);

%Taking the neighborhood around the window
idx = xx-windowSize:xx+windowSize;
idy = yy-windowSize:yy+windowSize;

%sChecking if out of bounds
idy(idy<1)=1;
idx(idx<1)=1;
idy(idy>size(filterResponse,2))=size(filterResponse,2);
idx(idx>size(filterResponse,1))=size(filterResponse,1);

%Assigning the small window (i.e. 11 by 11) around the peak as 0
filterResponse(idx,idy)=0;

%mean value and the standard deviation of the sidelobe
m = sum(filterResponse(:))/(numel(filterResponse)- (windowSize^2));
d=sqrt(size(filterResponse(:),1)*var(filterResponse(:))/(numel(filterResponse)- (windowSize^2)));
psr =(maxresponse - m)/d ;

end
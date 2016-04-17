
%  Exploiting the Circulant Structure of Tracking-by-detection with Kernels
%
%  Main script for tracking, with a gaussian kernel.
%
%  João F. Henriques, 2012
%  http://www.isr.uc.pt/~henriques/


%choose the path to the videos (you'll be able to choose one with the GUI)
base_path = './tiger2/';
video_path = './tiger2/imgs/';
close all; 

%parameters according to the paper
padding = 1;					%extra area surrounding the target
output_sigma_factor = 1/16;		%spatial bandwidth (proportional to target)
sigma = 0.2;					%gaussian kernel bandwidth
lambda = 1e-2;					%regularization
interp_factor = 0.075;			%linear interpolation factor for adaptation


clearvars A b C S Sp pos posvec

%notation: variables ending with f are in the frequency domain.

%ask the user for the video
%{
video_path = choose_video(base_path);
if isempty(video_path), return, end  %user cancelled
%}
[img_files, pos, target_sz, resize_image, ground_truth, video_path] = ...
	load_video_info(video_path);


%window size, taking padding into account
sz = floor(target_sz * (1 + padding));

%desired output (gaussian shaped), bandwidth proportional to target size
output_sigma = sqrt(prod(target_sz)) * output_sigma_factor;
[rs, cs] = ndgrid((1:sz(1)) - floor(sz(1)/2), (1:sz(2)) - floor(sz(2)/2));
y = exp(-0.5 / output_sigma^2 * (rs.^2 + cs.^2));
yf = fft2(y);

%store pre-computed cosine window
cos_window = hann(sz(1)) * hann(sz(2))';


time = 0;  %to calculate FPS
positions = zeros(numel(img_files), 2);  %to calculate precision
occluded = 0;
for frame = 1:numel(img_files),
	
    %% load image
	im = imread([video_path img_files{frame}]);
	if size(im,3) > 1,
		im = rgb2gray(im);
	end
	if resize_image,
		im = imresize(im, 0.5);
	end
	
	tic()
	
	%extract and pre-process subwindow
    if(occluded)
        x = get_subwindow(im, next_pos_vel, sz, cos_window);
    else
        x = get_subwindow(im, pos, sz, cos_window);
    end
    
    %% calculate the gaussian response
    if frame > 1,
        %calculate response of the classifier at all locations
        k = dense_gauss_kernel(sigma, x, z);
        response = real(ifft2(alphaf .* fft2(k)));   %(Eq. 9)
        psr = PSR(response);
        
        % determine if occluded
        occluded = (PSR(response) < 10);
        
        %target location is at the maximum response
        [row, col] = find(response == max(response(:)), 1);
        pos = pos - floor(sz/2) + [row, col];
        vel = pos - last_pos;
    end
    if(~occluded)
        %% adjust model
        %get subwindow at current estimated target position, to train classifer
        x = get_subwindow(im, pos, sz, cos_window);
        
        %Kernel Regularized Least-Squares, calculate alphas (in Fourier domain)
        k = dense_gauss_kernel(sigma, x);
        new_alphaf = yf ./ (fft2(k) + lambda);   %(Eq. 7)
        new_z = x;
        
        if frame == 1,  %first frame, train with a single image
            alphaf = new_alphaf;
            z = x;
            response = zeros(size(k));
            next_pos = pos;
            vel = [0 0];
            last_pos = pos;
            next_pos_vel = pos;
        else
            %subsequent frames, interpolate model only if not occluded
            alphaf = (1 - interp_factor) * alphaf + interp_factor * new_alphaf;
            z = (1 - interp_factor) * z + interp_factor * new_z;
        end
    end
    
    
    %% store the state values
    n = 8; d = 4;
    
    % create vector of states (positions), posvec, which is (1 x n x 2)
    if(frame <= 2*n)
        % initially just populate vector
        statevec(1,frame,:) = shiftdim([pos vel],-1);
    else
        % drop the first, add the last
        if(occluded)
            statevec = [statevec(1,2:end,:) shiftdim(next_state,-2)];
        else
            statevec = [statevec(1,2:end,:) shiftdim([pos vel],-1)];
        end
    end
    
    %% predict next position if occluded
    if(frame > 2*n)
        
        %         % smooth the data like a sly dog
        %         for m = 1:d
        %             statevec(:,:,m) = smooth(statevec(:,:,m),7);
        %         end
        %
        
        %% create square Hankel matrix
        Hxyd = zeros([n,n,d]);
        
        for ii = 1:n
            for jj = 1:n
                m = (ii-1)+ (jj-1) +1; % vector index
                Hxyd(ii,jj,:) = statevec(1,m,:); % Hankel my ankle
            end
        end
        
        % fold higher dimension vector into a block Hankel matrix
        Hxy = reshape(shiftdim(Hxyd,2),[],n,1);

        %% estimate complexity
        %         % build A matrix
        Axyd = Hxyd(:,1:end-1,:);
        %         bxy = Hxyd(:,end,:);
        %         Cxy = Hxyd(end,2:end,:);
        
        Axy = reshape(shiftdim(Axyd,2),[],n,1);
        
        S{frame} = svd(Axy);
        
        % decide complexity
        k = 3;
        %{
        [U, ~, V] = svd(A(:,:,m));
        S = svd(A(:,:,m));
        Sp = zeros(size(S));
        Sp(1:l) = S(1:l);
        Ap(:,:,m) = U*diag(Sp)*V';
        %}
        
        % decide memory
        nn = 2*n-k+1; % total available rows
        mm = min(k+2,nn); % number of rows to keep
        %% rebuild block hankel matrix to enforce minimum rank
        Hpd = zeros([mm,k,d]);
        for ii = nn-mm+1:nn % choose final mm rows
            for jj = 1:k
                m = (ii-1) + (jj-1) + 1; % vector index
                Hpd(ii-nn+mm,jj,:) = statevec(1,m,:); % Hankel my ankle
            end
        end

        %         % fold higher dimension vector into a block Hankel matrix
        Hp = reshape(shiftdim(Hpd,2),[],k,1);
        
        %% estimate next state
        % get sample vectors as a 2D vector matrix
        Ad = Hpd(:,1:end-1,:);
        bd = Hpd(:,end,:);
        Cd = Hpd(end,2:end,:);
        
        % fold vectors into 2 dimensional block Hankel matrix
        A = reshape(shiftdim(Ad,2),[],k-1,1);
        b = reshape(shiftdim(bd,2),[],1,1);
        C = reshape(shiftdim(Cd,2),[],k-1,1);
                
        % compute dynamic linear regressor coefficients, v
        v = A \ b;
        
        % predict next state (location)
        next_state = C*v;
        next_pos = next_state(1:2);
        next_pos_vel = (pos(:) + next_state(3:4)).';
        
    end
    
    %% save position and calculate FPS
	if(occluded)
        positions(frame,:) = next_pos;
    else
        positions(frame,:) = pos;
    end
    last_pos = pos;
	time = time + toc();
    
	%% visualization
	rect_position = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
    rect_2 = [next_pos_vel([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
	if frame == 1,  %first frame, create GUI
		f = figure('NumberTitle','off', 'Name',['Tracker - ' video_path]);
        im_handle = imshow(im, 'Border','tight', 'InitialMag',200);
		rect_handle = rectangle('Position',rect_position, 'EdgeColor','g');
        rect_handle2 = rectangle('Position',rect_2,'EdgeColor','r');
        ax = axis;
        [Ny, Nx] = size(response);
        [Xx, Yy] = meshgrid(linspace(rect_position(1),rect_position(1) + rect_position(3),Nx),...
                            linspace(rect_position(2),rect_position(2) + rect_position(4),Ny));
		hold on;
        p_handle = pcolor(Xx,Yy,128+128*response); 
        shading interp; 
        set(p_handle,'FaceAlpha','interp',...
            'AlphaData',response.^3,...
            'FaceColor','g');
        np_handle(1) = plot(next_pos(2),next_pos(1),'m+');
        np_handle(2) = plot(next_pos_vel(2),next_pos_vel(1),'cx');
        axis(ax);
        placeplot(7,f(1));
        f(2) = figure; placeplot(1,f(2)); axs2 = gca;
        f(3) = figure; placeplot(3,f(3)); axs3 = gca;
        
	else
		try  %subsequent frames, update GUI
			set(im_handle, 'CData', im)
			set(rect_handle, 'Position', rect_position)
            set(rect_handle2, 'Position', rect_2);
            [Xx, Yy] = meshgrid(linspace(rect_position(1),rect_position(1) + rect_position(3),Nx),...
                linspace(rect_position(2),rect_position(2) + rect_position(4),Ny));
            set(p_handle,'AlphaData',response.^3);
            set(p_handle,'Xdata',Xx,'Ydata',Yy);
            set(np_handle(1),'Xdata',next_pos(2),'Ydata',next_pos(1));
            set(np_handle(2),'Xdata',next_pos_vel(2),'Ydata',next_pos_vel(1));
            if(occluded)
                np_handle(1).Color = 'm';
                np_handle(2).Color = 'c';
            else
                np_handle(1).Color = 'r';
                np_handle(2).Color = 'b';
            end
            set(gcf,'name',(sprintf('PSR ~ %0.2f\tNext position at (%i, %i) ?',...
                PSR(response),round(next_pos(1)),round(next_pos(2)))));
            plot(axs2,statevec(:,:,1));
            plot(axs3,statevec(:,:,2));
            
        catch %user has closed the window
            return
		end
	end
	
	drawnow
%  	pause(0.15)  %uncomment to run slower
%     if(~mod(frame,10)); waitforbuttonpress; end
%     if(occluded); waitforbuttonpress; end
    if(frame == 104); waitforbuttonpress; end
end

if resize_image, positions = positions * 2; end

disp(['Frames-per-second: ' num2str(numel(img_files) / time)])

%show the precisions plot
show_precision(positions, ground_truth, video_path)


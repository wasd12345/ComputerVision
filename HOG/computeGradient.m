function [ angles, magnitudes ] = computeGradient( im )
%COMPUTEGRADIENT Given an image, computes the pixel gradients
% Arguments:
%          im - an image matrix read in by im2read (size H X W X C)
%
% Returns:
%          angles - (H-2) x (W-2) matrix of gradient angles
%          magnitudes - (H-2) x (W-2) matrix of gradient magnitudes

H = size(im,1);
W = size(im,2);
C = size(im,3);

%Empty container arrays:
stack3channel_x = zeros(H-2,W-2,C);
stack3channel_y = zeros(H-2,W-2,C);

%Iterate through each color channel:
for c=1:C
    
    %Get the various offsets of the pixels
    arr = im(:,:,c);
    from_right = circshift(arr,-1,2);
    from_left = circshift(arr,1,2);
    from_down = circshift(arr,-1,1);
    from_up = circshift(arr,1,1);    
    
    dx = from_left - from_right;
    dx = dx(2:end-1,2:end-1);
    
    dy = from_up - from_down;
    dy = dy(2:end-1,2:end-1);    

    %Insert into the container arrays:
    stack3channel_x(:,:,c) = dx;
    stack3channel_y(:,:,c) = dy;    

end


%Now get the channel-wise max for each pixel:
%[Alternative would be to use the color with the max gradient magnitude
%(instead of looking at the dx and dy components individually)]
dx_max = max(stack3channel_x,[],3);
dy_max = max(stack3channel_y,[],3);

%Magnitudes and angles
magnitudes = sqrt(dx_max.^2 + dy_max.^2);
angles = atand(dx_max./dy_max) + 90;

%To be safe, set any NANs or infs to 0
magnitudes(isnan(magnitudes)) = 0;
magnitudes(isinf(magnitudes)) = 0;
angles(isnan(angles)) = 0;
angles(isinf(angles)) = 0;    

showplot = false;
if showplot==true
    %Magnitude plot
    figure()
    imshow(mat2gray(magnitudes))
    colorbar()
    %Angle plot
    figure()
    imshow(mat2gray(angles))
    colormap('jet')
    colorbar()
    figure()
    title('Sorted Angles (degrees)')
    plot(sort(angles(:)))
end


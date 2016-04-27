function [F,d_ave] = GetF(correspondences, im1, im2, normalization)
% Use the 8 Point Algorithm (optionally normalized) on the 2 image
% correspondences to calculate the Fundamental Matrix.
% x'^T * F * x = 0 can be written as Af=0, where A is the Npoints x 9
% matrix of constraints, and f is the flattened 9 x 1 vector of F elements.
% Once F is known, compare l'=Fx to the point x' (which ideally is on l'),
% and also compare l=F^Tx' to x (which ideally is on l). Return the
% average distance between the points and their corresponding epipolar
% lines, and plot images with the points and epipolar lines overlaid.

% Format of correspondences is: size = Npoints x 4:
% [x1, y1, x2, y2]

format long

% Find number of point correspondences
N = size(correspondences,1);

% First check to make sure you have >= 8 points:
if N < 8
    error('Assumes >= 8 points to use the 8 point algorithm')
end


% NORMALIZATION:
%--------------------------------------------------------------------------
% Do a normalizing transformation
% x^ = Tx, x^' = Ux'
% where T and U are 3x3 transformation matrices that translate and
% scale the coordinates to be 0 mean, with average distance from origin 
% sqrt(2) for 2D image coordinates

if normalization == true
    
    % Mean of each coordinate vector
    mean_vector = mean(correspondences);

    % For 2D image1 coordinates, get average distance to centroid
    dx1 = correspondences(:,1) - mean_vector(1);
    dy1 = correspondences(:,2) - mean_vector(2);
    d1 = mean(sqrt(dx1.^2 + dy1.^2));
    K1 = sqrt(2)/d1;

    % For 2D image2 coordinates, get average distance to centroid
    dx2 = correspondences(:,3) - mean_vector(3);
    dy2 = correspondences(:,4) - mean_vector(4);
    d2 = mean(sqrt(dx2.^2 + dy2.^2));
    K2 = sqrt(2)/d2;
    
    % Make the transformation matrices
    T1 = [[K1, 0, -K1*mean_vector(1)]; [0, K1, -K1*mean_vector(2)]; [0, 0, 1]];
    T2 = [[K2, 0, -K2*mean_vector(3)]; [0, K2, -K2*mean_vector(4)]; [0, 0, 1]];

    % Apply the normalization transformation to image1 points
    trnsf1 = vertcat(correspondences(:,1).',correspondences(:,2).',ones(1,N));
    trnsf1 = (T1*trnsf1).'; %is column vectors of [x1, y1, 1]
    x1 = trnsf1(:,1);
    y1 = trnsf1(:,2);
    
    % Apply the normalization transformation to image2 points
    trnsf2 = vertcat(correspondences(:,3).',correspondences(:,4).',ones(1,N));
    trnsf2 = (T2*trnsf2).'; %is column vectors of [x2, y2, 1]
    x2 = trnsf2(:,1);
    y2 = trnsf2(:,2);
    
else
    % If not doing normalization (not advised)
    x1 = correspondences(:,1);
    y1 = correspondences(:,2);
    x2 = correspondences(:,3);
    y2 = correspondences(:,4);
    
end
%--------------------------------------------------------------------------



% BUILDING CONSTRAINT MATRIX:
%--------------------------------------------------------------------------
A = horzcat(x2.*x1, x2.*y1, x2, y2.*x1, y2.*y1, y2, x1, y1, ones(N,1));
%--------------------------------------------------------------------------



% SOLVING MATRIX EQUATION:
%--------------------------------------------------------------------------
[fU, fS, fV] = svd(A);
f = fV(:,end);
Ftemp = reshape(f,3,3).';
%--------------------------------------------------------------------------



% ENFORCING RANK 2 CONSTRAINT ON F:
%--------------------------------------------------------------------------
% Find SVD
[U, S, V] = svd(Ftemp);
% Set 3rd singular value to 0 to enforce rank 2
S(3,3) = 0;
% Recalculate F
F = U*S*V.';
%--------------------------------------------------------------------------



% DENORMALIZING:
%--------------------------------------------------------------------------
if normalization == true
    F = T2.'*F*T1;
end
%--------------------------------------------------------------------------



% ERROR ANALYSIS:
%--------------------------------------------------------------------------
x = vertcat(correspondences(:,1).',correspondences(:,2).',ones(1,N));
lp = F*x; %l'=Fx is the epipolar line corresponding to x
xp = vertcat(correspondences(:,3).',correspondences(:,4).',ones(1,N));
l = F.'*xp; %l=F^T x' is the epipolar line corresponding to x'

% Normalize homogeneous lines l and lp:
l = l./repmat(sqrt(l(1,:).^2+l(2,:).^2),3,1);
lp = lp./repmat(sqrt(lp(1,:).^2+lp(2,:).^2),3,1);

% Get perpendicular distances of points to corresponding epipolar lines:
d1_ave = mean(abs(dot(x,l)));
d2_ave = mean(abs(dot(xp,lp)));
% Average over both images (which of course have same number of points
% since they are using the same 3D points) to get single average error for
% the data set:
d_ave = mean([d1_ave,d2_ave]);
%--------------------------------------------------------------------------



% PLOTTING POINTS AND EPIPOLAR LINES:
%--------------------------------------------------------------------------
%ax+by+c=0 -> y=(-c-ax)/b
%So endpoints: [x-d,(-c-a(x-d))/b] and [x+d,(-c-a(x+d))/b]

%size of increment d:
d = 20;

%epipolar lines l vs. x
a = l(1,:);
b = l(2,:);
c = l(3,:);
x_minus_d = x(1,:) - repmat(d,1,N);
x_plus_d = x(1,:) + repmat(d,1,N);
y_minus = (-c - a.*x_minus_d)./b;
y_plus = (-c - a.*x_plus_d)./b;

%epipolar lines l' vs. x' (p suffix meaning prime)
ap = lp(1,:);
bp = lp(2,:);
cp = lp(3,:);
xp_minus_d = xp(1,:) - repmat(d,1,N);
xp_plus_d = xp(1,:) + repmat(d,1,N);
yp_minus = (-cp - ap.*xp_minus_d)./bp;
yp_plus = (-cp - ap.*xp_plus_d)./bp;


% Plot images, overlaid with points (blue) and epipolar lines (red)
figure
% LEFT: im1 w/ l and x,
subplot(1,2,1), imshow(im1); hold on; plot(x(1,:), x(2,:), '.', 'markers',12);
hold on; plot(vertcat(x_minus_d,x_plus_d), vertcat(y_minus,y_plus), 'r');
% RIGHT: im2 w/ l' and x'
subplot(1,2,2), imshow(im2); hold on; plot(xp(1,:), xp(2,:), '.', 'markers',12);
hold on; plot(vertcat(xp_minus_d,xp_plus_d), vertcat(yp_minus,yp_plus), 'r');
%--------------------------------------------------------------------------
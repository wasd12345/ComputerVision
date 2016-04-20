function [P_affine, RMSerror, projected] = AffineCamera(coordinates)
% Calculate Affine Camera Matrix from correspondences

% Calculate affine camera matrix, P_affine,
% and RMS error of N image point locations <-> N calculated locations.
% Inputs: coordinates is Nx5 array of image and world point coordinates.
% Format is: [X_world, Y_world, Z_world, x_image, y_image]
% Basic idea is to solve for p from the equation: Ap = b

% the output "projected" is the x,y image coordinates predicted by using
% the affine camera matrix on the world coordinates of the calibration
% grid. It's useful for error analysis.

% Find number of point correspondences
N = size(coordinates,1);


% NORMALIZATION:
%--------------------------------------------------------------------------
%Do a normalizing transformation
% x' = Tx, X' = UX
% where T is 3x3, U is 4x4, are transformation matrices that translate and
% scale the coordinates to be 0 mean, with average distance from origin 
% sqrt(2) for 2D image coordinates and sqrt(3) for 3D world coordinates
mean_vector = mean(coordinates); %mean of each coordinate vector

% For #D world coordinates, get average distance to centroid
dX = coordinates(:,1) - mean_vector(1);
dY = coordinates(:,2) - mean_vector(2);
dZ = coordinates(:,3) - mean_vector(3);
d_ave__3D = sum(sqrt(dX.^2 + dY.^2 + dZ.^2))/N;
K_3D = sqrt(3)/d_ave__3D;

% For 2D image coordinates, get average distance to centroid
dx = coordinates(:,4) - mean_vector(4);
dy = coordinates(:,5) - mean_vector(5);
d_ave__2D = sum(sqrt(dx.^2 + dy.^2))/N;
K_2D = sqrt(2)/d_ave__2D;

% Make the transformation matrices
T = [[K_2D, 0, -K_2D*mean_vector(4)]; [0, K_2D, -K_2D*mean_vector(5)]; [0, 0, 1]];
U = [[K_3D, 0, 0, -K_3D*mean_vector(1)]; [0, K_3D, 0, -K_3D*mean_vector(2)]; [0, 0, K_3D, -K_3D*mean_vector(3)]; [0, 0, 0, 1]];

% Apply the normalization transformation to world coordinates
trnsf3D = vertcat(coordinates(:,1).',coordinates(:,2).',coordinates(:,3).',ones(1,N));
trnsf3D = (U*trnsf3D).'; %is column vectors of [X, Y, Z, 1]

% Apply the normalization transformation to image coordinates
trnsf2D = vertcat(coordinates(:,4).',coordinates(:,5).',ones(1,N));
trnsf2D = (T*trnsf2D).'; %is column vectors of [x, y, 1]
%--------------------------------------------------------------------------




% BUILDING MATRICES AND SOLVING SYSTEM OF EQUATIONS:
%--------------------------------------------------------------------------
% Build a 2N x 8 matrix representing constraints from each of the N points
top_half = horzcat(trnsf3D(:,1:3), repmat([1, 0, 0, 0, 0], N, 1));
bottom_half = horzcat(zeros(N,4), trnsf3D);
A = vertcat(top_half,bottom_half);

% Build 2N x 1 vector of transformed x, then y image points
b = vertcat(trnsf2D(:,1),trnsf2D(:,2));

% Solve Ap = b for p vector (camera matrix parameters)
p8 = pinv(A)*b;
P = vertcat(p8(1:4).',p8(5:8).',[0,0,0,1]);

% Since normalization was used before, now denormalize it
P_affine = T\(P*U); %T is invertible.
%--------------------------------------------------------------------------





% ERROR ANALYSIS:
%--------------------------------------------------------------------------
% Calculate RMS error of N image point locations <-> N calculated locations
% Get the projected coordinates using the affine camera matrix
projected = P_affine*vertcat(coordinates(:,1:3).',ones(1,N));
%Since affine camera, last row of P_affine is [0,0,0,1], so to go from
% homogeneous image coordinates -> Euclidean, last element is 1, so don't
% need to divide: can just drop last element:
projected = projected(1:2,:);
% projected is a 2 x N array w/ format: [x_projected; y_projected]
% the actual x and y to compare to arethe measured image coordinates:
x_diff = coordinates(:,4).' - projected(1,:);
y_diff = coordinates(:,5).' - projected(2,:);
x_squares = x_diff.^2;
y_squares = y_diff.^2;
x_sum_squares = sum(x_squares(:));
y_sum_squares = sum(y_squares(:));
% RMSerror
RMSerror = sqrt((x_sum_squares + y_sum_squares)/N);
%--------------------------------------------------------------------------
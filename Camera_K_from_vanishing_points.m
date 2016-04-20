% Find focal length and camera center using vanishing points.

% OVERVIEW
% -------------------------------------------------------------------------
% Overview: following algorithm 8.2 on HZ pg. 225, and HZ pg. 226.
% Omega, the image of the absolute conic (IAC), is a 3x3 matrix. Any 3x3
% symetric matrix has 6 degrees of freedom. However, since we are told that
% the camera has 0 skew, and square pixels (unit aspect ratio), there are 
% now only 4 elements of omega to solve for. Using 3 orthogonal vanishing 
% points, we can solve for the elements of omega.
% Once we have omega, we can use a simple relationship between omega and
% the camera calibration matrix, K, to determine K, and then look at
% camera parameters.
% -------------------------------------------------------------------------

%change formatting options so can see values better
format long





% LOAD VANISHING POINTS
% -------------------------------------------------------------------------
% Load set of 3 orthogonal vanishing pts. obtained from getVanishingPoint.m
% Is [[x1, y1, 1]; [x2, y2, 1]; [x3, y3, 1]], is a 3x3 array of
% homogeneous image coordinates of vanishing points:
% top left from floor tiles; top right from floor; bottom from shelf
% is [[ -707.2465, -120.2651, 1]; [5.1073e+03, -0.0755e+03, 1]; 
% [1.2810e+03, 3.9331e+03, 1]]
v = load('orthogonal_vanishing_points.mat');
v = v.orthogonal_vanishing_points;
% One of the vanishing points (vertical) was bad, so use a better one:
%image1 cardboard box: use left edge and edge center to find vertical VP:
v(3,:) = [1.153141892379285e03, 5.698086563406738e03, 1];
display 'Vanishing Points used:'
v
% -------------------------------------------------------------------------







% CREATE SYSTEM OF EQUATIONS
% -------------------------------------------------------------------------
% Make stack of 3 constraint equations from orthogonal vanishing points
% Derivation is provided on separate page
% A is a 3x4 matrix, imposing 3 constraints on omega
A11 = v(1,1)*v(2,1) + v(1,2)*v(2,2);
A21 = v(1,1)*v(3,1) + v(1,2)*v(3,2);
A31 = v(2,1)*v(3,1) + v(2,2)*v(3,2);
A12 = v(1,1) + v(2,1);
A22 = v(1,1) + v(3,1);
A32 = v(2,1) + v(3,1);
A13 = v(1,2) + v(2,2);
A23 = v(1,2) + v(3,2);
A33 = v(2,2) + v(3,2);
A = [[A11,A12,A13,1]; [A21,A22,A23,1]; [A31,A32,A33,1]];
% -------------------------------------------------------------------------





% SOLVE SYSTEM OF EQUATIONS
% -------------------------------------------------------------------------
% Solve matrix equation A*w = 0 for w (4x1 vector of the elements of omega)
% by getting right singular vector corresponding to smallest singular value
[U, S, V] = svd(A);
w = V(:,end);
% -------------------------------------------------------------------------





% MAKE OMEGA (IAC) MATRIX, AND K (CAMERA CALIBRATION MATRIX)
% -------------------------------------------------------------------------
% Construct omega matrix from the elements w:
omega = [[w(1),0,w(2)]; [0,w(1),w(3)]; [w(2),w(3),w(4)]]
% check omega is positive definite:
eigs = eig(omega);

% Cholesky factorization to get K:
%K = chol(inv_omega);
chol_omega = chol(omega);
K = inv(chol_omega);

% Scale K to have K(3,3) = 1 to follow convention
K = K/K(3,3)

% Extract the focal length from K:
% Since unit aspect ratio, should have identical K(1,1) and K(2,2),
% so can use either, or average both elements just in case:
f = (K(1,1) + K(2,2)) / 2

% Calculate the camera center:
% The principal offset is (K(1,3), K(2,3), and is the projection of camera
% center on the image plane. The camera center has the same x, y
% coordinates but has a z coordinate that is offset by the focal length. In
% the commonly used right hand coordinate system where the image plane is
% in front of the pinhole (same side as incoming light) and with z
% increasing (+) from image plane to pinhole, we have:
principal_point = [K(1,3); K(2,3)]
disp 'Camera Center'
C = [K(1,3); K(2,3); f]
% -------------------------------------------------------------------------
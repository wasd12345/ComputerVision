% Camera is rotated (but not translated) to get from image1 to image2.
% K matrix remains the same (intrinsics unchanged).
% Use vanishing points to find rotation matrix.


%--------------------------------------------------------------------------
% Overview: vanishing points are images of points at infinity.
% So vanishing points are unchanged by camera translation, but rotation does 
% change the vanishing points. Vanishing points can be related to directions
% through the matrix K, and each direction d_i from image 1 is related to 
% a direction d_i' in image2 by a rotation:
% d_i' = R * d_i

% Each direction has 3 components, so each direction provides 3 equations
% for elements of R. R is a 3x3 rotation matrix so has 9 elements.

% So by stacking 3 direction relationships together (by using the 
% corresponding vanishing points from image1 and image2), you will have 9
% equations for elements of R and can estimate the rotation matrix by 
% inverting a matrix. You could use more than 3 directions, and create an 
% over-determiend system and use least squares, but the given set of images
% acquireed in the lab of the cardboard box basically only have 3 strong
% othogonal directions, (the scene is basically a 3-point perspective),so 
% I chose to just use the best vanishing points I could to get reliable 
% direction vectors and solve the system exactly.

% di = K^-1 * vi / ||K^-1 * vi||, and similar for di' <-> vi' relationship
%--------------------------------------------------------------------------



% Change formatting options so can see values better
% format long

% Define K, the matrix of camera intrinsics, given as:
K = [[2448, 0, 1253]; [0, 2438, 986]; [0, 0, 1]]






% LOAD DATA
%--------------------------------------------------------------------------
% Concatenate the vanishing points from image1:
v_im1 = load('Q3C_im1_points.mat');
display 'image1 vanishing points: each column is: (x,y,1)^T'
v_im1 = v_im1.points %is 3x3 array of: each column is a vanishing point
% make d_im1, the 3x4 array where each column is a (unit) direction in im1
%display 'image1 unit direction vectors: each column is: (x,y,z)^T'
d_im1 = inv(K)*v_im1 ./ repmat(sqrt(sum((inv(K)*v_im1).^2)),3,1); %(MATLAB default sum is columnwise)

% Do the same for the image2 vanishing points to get directions in image2:
%v_im2 = load('Q3C_im2_points.mat');
%v_im2 = v_im2.points2 %is 3x3 array of: each column is a vanishing point
% Using a better set of points in image2:
v_im2 = load('Q3C_im2_points3.mat');
display 'image2 vanishing points: each column is: (x,y,1)^T'
v_im2 = v_im2.points3 %is 3x3 array of: each column is a vanishing point
%display 'image2 unit direction vectors: each column is: (x,y,z)^T'
d_im2 = inv(K)*v_im2 ./ repmat(sqrt(sum((inv(K)*v_im2).^2)),3,1);
%--------------------------------------------------------------------------




% CREATE MATRIX EQUATIONS
%--------------------------------------------------------------------------
% Create 9x1 solution vector, made of all 3 components of 3 di' vectors
b = d_im2(:);

% Create constraint matrix from components of di vectors
A1 = horzcat(d_im1(:,1).',zeros(1,6));
A2 = circshift(A1,3,2);
A3 = circshift(A2,3,2);

A4 = horzcat(d_im1(:,2).',zeros(1,6));
A5 = circshift(A4,3,2);
A6 = circshift(A5,3,2);

A7 = horzcat(d_im1(:,3).',zeros(1,6));
A8 = circshift(A7,3,2);
A9 = circshift(A8,3,2);

A = vertcat(A1,A2,A3,A4,A5,A6,A7,A8,A9);
%--------------------------------------------------------------------------



% SOLVE SYSTEM OF EQUATIONS
%--------------------------------------------------------------------------
% Look at singular values to see how singular the matrix is
[U, S, V] = svd(A);

% Solve the system exactly for r (elements of rotation matrix R):
%A is an invertible matrix, so pseudo-inverse gives same result as inv
%r = A \ b;
r = pinv(A)*b; 

% reshape to get matrix R
R = reshape(r,3,3).'
%--------------------------------------------------------------------------




% BASIC VERIFICATION OF RESULTS
%--------------------------------------------------------------------------
% Basic checks of the rotation matrix:
% ideally, it'san orthogonal matrix (orthonormal column vectors), with
% would have determinant +1 (+ since non-reflective), 
% and R*R^T = I (Identity matrix 3x3)
detR = det(R)
RRT = R*R.'
display 'Dot products should be 0 (orthogonal columns)'
u1u2 = dot(R(:,1),R(:,2))
u1u3 = dot(R(:,1),R(:,3))
u2u3 = dot(R(:,2),R(:,3))
display 'Norms should be 1 (orthonormal columns)'
norm1 = norm(R(:,1))
norm2 = norm(R(:,2))
norm3 = norm(R(:,3))
%--------------------------------------------------------------------------
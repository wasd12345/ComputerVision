function TomasiKanadeFactorization(usesubset)
% Use the Tomasi Kanade Factorization method to
% solve the affine structure from motion problem


%Read in image files
[x1, x2, pts3d] = readTextFiles(usesubset);
x1=x1(1:2,:);
x2=x2(1:2,:);

%Whether or not to use subset of points:
usesubset
Npoints = size(x1,2);

% CENTERING THE DATA:
%--------------------------------------------------------------------------
%Center the image coordinates for image1 and image2:
x1 = x1 - repmat(mean(x1,2),1,Npoints);
x2 = x2 - repmat(mean(x2,2),1,Npoints);
%--------------------------------------------------------------------------


% CREATE MEASUREMENT MATRIX D:
%--------------------------------------------------------------------------
% D is 2m x n, where m=Ncameras, n=Npoints
D = vertcat(x1,x2);
%--------------------------------------------------------------------------


% FACTORIZE THE MATRIX:
%--------------------------------------------------------------------------
[U,S,V] = svd(D);
disp 'Q3 c)and e): 4 Nonzero singular values:'
S4 = [S(1,1), S(2,2), S(3,3), S(4,4)]
%Have rank(D) > 3 since noisy measurements and affine approximation is 
%inexact, so enforce rank 3 via factorization:
%Structure and motion matrices:
STRUCTURE = S(1:3,1:3)*V(:,1:3).';
MOTION = U(:,1:3);
%Check reasonable approximation:
diff = D - MOTION*STRUCTURE;
%--------------------------------------------------------------------------


% PLOT THE RESULTING 3D POINTS:
%--------------------------------------------------------------------------
figure()
plot3(STRUCTURE(1,:),STRUCTURE(2,:),STRUCTURE(3,:),'ro')
axis equal
title('Calculated World Points')

figure()
plot3(pts3d(1,:),pts3d(2,:),pts3d(3,:),'bo')
axis equal
title('Ground Truth World Points')
%--------------------------------------------------------------------------


end
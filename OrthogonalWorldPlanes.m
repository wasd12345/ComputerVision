% Change formatting options so can see values better
format long

% Define K, the matrix of camera intrinsics, given as:
K = [[2448, 0, 1253]; [0, 2438, 986]; [0, 0, 1]]

% Load set of 4 vanishing points
disp 'vanishing points'
v = load('3b_points.mat');
v = v.points
% row 1 is vanishing point of box at bottom
% row 2 is vanishing point of box to left
% row 3 is the vanishing point of the ground to left
% row 4 is the vanishing point of the ground to right


% Verify that the ground plane is orthogonal to the plane on the cardboard 
% box containing the letter "A"
% Get the normal vector for plane A on the box:
% identify the horizon line for plane A by connecting 2 vanishing points.
% those 2 vanishing points correspond to the 2 vertical parallel lines on 
% plane A, and the 2 horizontal lines of plane A.


% Plane A vertical vanishing point:
v1 = v(1,:).';
% plane A horizontal vanishing point:
v2 = v(2,:).';
% Get vanishing line
l1 = cross(v1,v2);

% Do same for ground plane
% ground plane left vanishing point:
v3 = v(3,:).';
% ground plane right vanishing point:
v4 = v(4,:).';
% Get vanishing line
l2 = cross(v3,v4);


% Get normal lines
n1 = K.' * l1;
n2 = K.' * l2;
% normalize to get unit vectors:
n1 = n1 / norm(n1);
n2 = n2 / norm(n2);


% Now verify that the planes are orthogonal:
% the dot product of 2 orthogonal vectors is 0:
deg2rad = pi/180;
dotprod = dot(n1,n2);
% dot(a,b) = ||a|| ||b|| cos(theta), and we have a and b unit vectors
theta = acos(dotprod)/deg2rad
disp 'Is ~ 90 degrees, meaning orthogonal'





% Alternative method as validation
% Using direction vectors from vanishing points
disp 'Alternative method as validation: use v = K*d -> d = K^-1 v'

% the 2 (normalized) direction vectors for vanishing points of plane A:
d1 = inv(K)*v1 / norm(inv(K)*v1); % direction of vertical vanishing point
d2 = inv(K)*v2 / norm(inv(K)*v2); % direction of plane A left vanishing pt.
d3 = inv(K)*v3 / norm(inv(K)*v3); % direction of ground left vanishing pt.
d4 = inv(K)*v4 / norm(inv(K)*v4); % direction of ground right vanishing pt.

% Get surface normal vectors
% sign (reflection) determines if facing in/out from surface
% but not important if we just want to check angle btwn. 2 surfaces
n1 = cross(d1,d2); %surface normal vector for plane A
n2 = cross(d3,d4); %surface normal vector for ground plane
% Normalize the vectors
n1 = n1/norm(n1); % surface normal unit vector for plane A
n2 = n2/norm(n2); % surface normal unit vector for ground plane

% Cross product of 2 perpendicular UNIT vectors is 1.
% ||cross(a,b)|| = ||a|| ||b|| sin(theta), and we have a and b unit vectors
c = cross(n1,n2);
cross_norm = norm(c); % So should be close to 1:
theta = asin(cross_norm)/deg2rad
disp 'Is ~ 90 degrees, meaning orthogonal'

% dot(a,b) = ||a|| ||b|| cos(theta), and we have a and b unit vectors
theta = acos(dotprod)/deg2rad
function [x1, x2, pts3d] = readTextFiles(usesubset)


% These are the points to sample when use_subset is set to true.
subset = [2 4 5 6 7 8 9 17 25 26 27 28 29 30 32 33 34 37];


% Read in the points for the first image.
fileid = fopen('set1/pt_2D_1.txt','r');
PointNo = str2num(fgets(fileid));
temp = fscanf(fileid,'%f');
Points1 = zeros(PointNo, 2);
for i = 1 : PointNo
    Points1(i,2) = temp(i*2-1);
    Points1(i,1) = temp(i*2);
end
fclose(fileid);

% Read in the points for the second image.
fileid = fopen('set1/pt_2D_2.txt','r');
PointNo = str2num(fgets(fileid));
temp = fscanf(fileid,'%f');
Points2 = zeros(PointNo, 2);
for i = 1 : PointNo
    Points2(i,2) = temp(i*2-1);
    Points2(i,1) = temp(i*2);
end
fclose(fileid);

if (usesubset)
    Points1 = Points1(subset,:);
    Points2 = Points2(subset,:);
end

figure; subplot(1,2,1), imshow('set1/image1.jpg'); hold on; plot(Points1(:,1), Points1(:,2), '.')
subplot(1,2,2), imshow('set1/image2.jpg'); hold on; plot(Points2(:,1), Points2(:,2), '.')

x1 = Points1';
x2 = Points2';

% Make into homogenous coordinates.
n = size( x1, 2 );
x1 = [ x1 ; ones(1,n) ];
x2 = [ x2 ; ones(1,n) ];

%Load ground truth
pts3d=importdata('set1/pt_3D.txt',' ',1);
pts3d=pts3d.data';

if (usesubset)
    pts3d = pts3d(:, subset);
end


end

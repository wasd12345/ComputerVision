Histogram of Oriented Gradients (HOG) image detector/descriptor

To run:
2 lines:

e.g.:


hogfeats = computeHOGFeatures(im2single(imread('car.jpg')), 8, 2, 9);
show_hog(hogfeats);
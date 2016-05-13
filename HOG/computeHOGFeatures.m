function [ features ] = computeHOGFeatures( im, cell_size, Ncells_per_block_dim, nbins )
%COMPUTEHOGFEATURES Computes the histogram of gradients features
% Arguments:
%          im - the image matrix
%          cell_size - each cell will be of size (cell_size, cell_size)
%                       pixels
%          Ncells_per_block - each block will be of size (Ncells_per_block_dim, Ncells_per_block_dim)
%                       cells
%          nbins - number of histogram bins
% Returns:
%          features - the hog features of the image (H_blocks x W_blocks x Ncells_per_block_dim*Ncells_per_block_dim*nbins)

% 0) Normalize the image globally (gamma and color)
% 0) Smooth the gradient image by convolving with Gaussian

% 1) Compute the gradient for the image
[grad_angle, grad_magnitude] = computeGradient(im);
H = size(grad_magnitude,1); %268 for car
W = size(grad_magnitude,2); %498 for car


% 2) Iterate through each block
block_pixel_width = cell_size*Ncells_per_block_dim;

%Feature vector container array
H_blocks = floor(H/(block_pixel_width/2)) - 1;
W_blocks = floor(W/(block_pixel_width/2)) - 1;

h=1; %counters for indices. Initially implemented features as 1D vector,
% so didn't need this, but realized later features was valency 3 tensor, so
% did this ad hoc.

features = zeros(H_blocks, W_blocks, Ncells_per_block_dim*Ncells_per_block_dim*nbins);
for yb=1:(block_pixel_width/2):(H-block_pixel_width)
    w=1;%counter
    for xb=1:(block_pixel_width/2):(W-block_pixel_width)
        %Technically don't need block-level data, only need the cell level,
        %but sttructuring code this way makes less error prone and easier to debug
        yb_end = yb+block_pixel_width-1;
        xb_end = xb+block_pixel_width-1;
        block_mags = grad_magnitude(yb:yb_end,xb:xb_end);
        block_angs = grad_angle(yb:yb_end,xb:xb_end);     
        
        %Iterate through the smaller cells within the block
        block_features = [];
        for yc=1:cell_size:block_pixel_width
            for xc=1:cell_size:block_pixel_width
                %Values of gradient magnitudes and angles in the cell
                %xc=xc
                %yc=yc
                yc_end = yc+cell_size-1;
                xc_end = xc+cell_size-1;
                cell_mags = block_mags(yc:yc_end,xc:xc_end);
                cell_mags = cell_mags(:);
                cell_angs = block_angs(yc:yc_end,xc:xc_end);
                cell_angs = cell_angs(:);

                %Histogram bins:
                bin_vals = zeros(1,nbins);
                %Input angles will be in range 0-180 degrees
                dtheta = 180./nbins;
                bin_boundaries = 0:dtheta:180;
                bin_centers = bin_boundaries(1:end-1)+dtheta/2;
                
                %Iterate through all of the pixels in the cell
                % Then share across 2 neighboring bins by:
                % |---x-----------|
                % |-a-|-----b-----|
                % give b/(a+b) to left, a/(a+b) to right
                % If is first or last bin need to wrap around.
                for pix=1:length(cell_mags)
                    ang = cell_angs(pix);
                    mag = cell_mags(pix);
                    %If angle is less than leftmost bin center:
                    if ang < bin_centers(1)
                        a = ang + (180-bin_centers(end));
                        b = bin_centers(1) - ang;
                        bin_vals(1) = bin_vals(1) + mag*(a/(a+b));
                        bin_vals(end) = bin_vals(end) + mag*(b/(a+b));
                    %If angle is greater than rightmost bin center:
                    elseif ang > bin_centers(end)
                        a = ang - bin_centers(end);
                        b = (180-ang) + bin_centers(1);
                        bin_vals(1) = bin_vals(1) + mag*(a/(a+b));
                        bin_vals(end) = bin_vals(end) + mag*(b/(a+b));
                    %If angle is within min and max bin centers:
                    else
                        for ind=1:nbins-1
                            c1 = bin_centers(ind);
                            c2 = bin_centers(ind+1);
                            if ((ang>=c1)&&(ang<c2))
                                a = ang-c1;
                                b = c2-ang;
                                bin_vals(ind) = bin_vals(ind) + mag*(b/(a+b));
                                bin_vals(ind+1) = bin_vals(ind+1) + mag*(a/(a+b));
                            end
                        end
                    end
                end
                block_features = horzcat(block_features,bin_vals);
            end
        end
        
        % Normalizae the block/cell (big structure)
        % Only do this if nonzero (important for Q3 C/D w. SVW weights)
        if sqrt(sum(block_features.^2)) ~= 0
            block_features = block_features/sqrt(sum(block_features.^2));
        end
        % Append to whole image feature tensor
        block_features = reshape(block_features,1,1,length(block_features));
        features(h,w,:) = block_features;
        w=w+1;%Increment counter
    end  
    h=h+1;%Increment counter
end

end

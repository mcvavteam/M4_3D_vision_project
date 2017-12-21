
function [image_transformed,transformed_corners] = apply_H (image, H)
% Apply the transformation H to the given image
% - image: nxmxc array, it can be a grayscale image or RGB image
% - H: transformation to apply to the given image

    % Obtain the image size
    [rows, cols, chan] = size(image);

    % Image corners in homogeneous coordinates.
    corners = zeros(3,4);
    corners(:,1) = [1 1 1]';
    corners(:,2) = [cols 1 1]';
    corners(:,3) = [1 rows 1]';
    corners(:,4) = [cols rows 1]';


    % Transform the image corners 
    transformed_corners = H*corners;

    % Make transformed_corners homogeneous
    for i=1:4
        transformed_corners(:,i) = transformed_corners(:,i)/transformed_corners(3,i);
    end

    % Compute the new extrenal transformed corner coordinates      
    xmin = round(min(transformed_corners(1,:)));
    xmax = round(max(transformed_corners(1,:)));
    ymin = round(min(transformed_corners(2,:)));
    ymax = round(max(transformed_corners(2,:)));

    % Create a grid with the new image size
    [X,Y] = meshgrid(xmin:xmax, ymin:ymax);
    Ncols = xmax - xmin + 1;
    Nrows = ymax - ymin + 1;
    W = ones(Nrows, Ncols);
    XYW = [X(:) Y(:) W(:)]';

    %Transform the image. p'= H*p --> p = inv(H)*p'
    transformed_XYW = H\XYW;

    transformed_X = reshape(transformed_XYW(1,:), Nrows, Ncols);
    transformed_Y = reshape(transformed_XYW(2,:), Nrows, Ncols);
    transformed_W = reshape(transformed_XYW(3,:), Nrows, Ncols);

    %Make transformed_X and transformed_Y homogeneous
    transformed_X = transformed_X./transformed_W;
    transformed_Y = transformed_Y./transformed_W;

    % Output to Input: Interpolate transformed points
    image_transformed = zeros(Nrows, Ncols, chan);

    for c=1:chan    
        image_transformed(:,:,c) = interp2(double(image(:,:,c)), ...
                                   transformed_X, transformed_Y, ...
                                   'linear', 0);
    end
end
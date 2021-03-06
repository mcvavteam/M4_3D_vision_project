%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 3: The geometry of two views 
% (application: photo-sequencing)

addpath('sift'); % ToDo: change 'sift' to the correct path where you have the sift functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Compute the fundamental matrix

% Two camera matrices for testing purposes
P1 = eye(3,4);
c = cosd(15); s = sind(15);
R = [c -s 0; s c 0; 0 0 1];
t = [.3 0.1 0.2]';
P2 = [R t];
n = 8;
X = [rand(3,n); ones(1,n)] + [zeros(2,n); 3 * ones(1,n); zeros(1,n)];
x1_test = P1 * X;
x2_test = P2 * X;

% Estimated fundamental matrix
% ToDo: create the following function that estimates F using the normalised 8 point algorithm
F_es = fundamental_matrix(x1_test, x2_test);

% Real fundamental matrix
T = [ 0  -t(3) t(2);
     t(3)  0  -t(1);
    -t(2) t(1)  0];

F_gt = T*R; % ToDo: write the expression of the real fundamental matrix for P1 and P2

% Evaluation: these two matrices should be very similar
F_gt / norm(F_gt)
F_es / norm(F_es)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Robustly fit fundamental matrix

% Read images
im1rgb = imread('Data/0000_s.png');
im2rgb = imread('Data/0001_s.png');
im1 = sum(double(im1rgb), 3) / 3 / 255;
im2 = sum(double(im2rgb), 3) / 3 / 255;

% show images
figure;
subplot(1,2,1); imshow(im1rgb); axis image; title('Image 1');
subplot(1,2,2); imshow(im2rgb); axis image; title('Image 2');

%% Compute SIFT keypoints

% (make sure that the sift folder provided in lab2 is on the path)

[points_1, desc_1] = sift(im1, 'Threshold', 0.01);
[points_2, desc_2] = sift(im2, 'Threshold', 0.01);

%% Match SIFT keypoints between a and b
matches = siftmatch(desc_1, desc_2);
figure;
plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches, 'Stacking', 'v');

% p1 and p2 contain the homogeneous coordinates of the matches
p1 = [points_1(1:2, matches(1,:)); ones(1, length(matches))];
p2 = [points_2(1:2, matches(2,:)); ones(1, length(matches))];

%% ToDo: create this function (you can use as a basis 'ransac_homography_adaptive_loop.m')
[F_12, inliers] = ransac_fundamental_matrix(p1, p2, 2.0, 1000); 
% F = fundamental_matrix(p1, p2); 
% inliers = 1:size(matches,2);

% show inliers
figure;
plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches(:,inliers), 'Stacking', 'v');
title('Inliers');

vgg_gui_F(im1rgb, im2rgb, F_12');

%% Plot some epipolar lines

l2 = F_12*p1 % epipolar lines in image 2 % ToDo
l1 = F_12'*p2 % epipolar lines in image 1 % ToDo

% choose three random indices
inliersSz = numel(inliers);
m1 = inliers(unidrnd(inliersSz));
m2 = inliers(unidrnd(inliersSz));
m3 = inliers(unidrnd(inliersSz));
% m1 = inliers(10);
% m2 = inliers(20);
% m3 = inliers(30);

% image 1 (plot the three points and their corresponding epipolar lines)
figure;
imshow(im1rgb);
hold on;
plot(p1(1, m1), p1(2, m1), '+g');
plot_homog_line(l1(:, m1));

plot(p1(1, m2), p1(2, m2), '+g');
plot_homog_line(l1(:, m2));

plot(p1(1, m3), p1(2, m3), '+g');
plot_homog_line(l1(:, m3));

% image 2 (plot the three points and their corresponding epipolar lines)
figure;
imshow(im2rgb);
hold on;
plot(p2(1, m1), p2(2, m1), '+g');
plot_homog_line(l2(:, m1));

plot(p2(1, m2), p2(2, m2), '+g');
plot_homog_line(l2(:, m2));

plot(p2(1, m3), p2(2, m3), '+g');
plot_homog_line(l2(:, m3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Photo-sequencing with aerial images

% In this part we will compute a simplified version of the algorithm
% explained in the Photo-sequencing paper. 
% Since we do not have two images
% taken from roughly the same viewpoint at two different time instants we
% will manually pick a dynamic point corresponding to a point in a van 
% (identified by index 'idx_car_I1') and the projection of its 3D trajectory 
% in the reference image. Then we will compute the projection (to the reference image) 
% of three points on this 3D trajectory at three different time instants 
% (corresponding to the time when the three other provided images where taken). 

clear all;

% Read images
im1rgb = imread('Data/frame_00000.tif');
im2rgb = imread('Data/frame_00001.tif');
im3rgb = imread('Data/frame_00002.tif');
im4rgb = imread('Data/frame_00003.tif');

im1 = sum(double(im1rgb), 3) / 3 / 255;
im2 = sum(double(im2rgb), 3) / 3 / 255;
im3 = sum(double(im3rgb), 3) / 3 / 255;
im4 = sum(double(im4rgb), 3) / 3 / 255;

% show images
figure;
subplot(2,2,1); imshow(im1rgb); axis image; title('Image 1');
subplot(2,2,2); imshow(im2rgb); axis image; title('Image 2');
subplot(2,2,3); imshow(im3rgb); axis image; title('Image 3');
subplot(2,2,4); imshow(im4rgb); axis image; title('Image 4');

% Compute SIFT keypoints
[points_1, desc_1] = sift(im1, 'Threshold', 0.015); % Do not change this threshold!
[points_2, desc_2] = sift(im2, 'Threshold', 0.015);
[points_3, desc_3] = sift(im3, 'Threshold', 0.015);
[points_4, desc_4] = sift(im4, 'Threshold', 0.015);

%% ToDo:

% Take image im1 as reference image (image 1) and compute the fundamental 
% matrices needed for computing the trajectory of point idx_car_I1
% (use the SIFT keypoints previously computed)

matches_12 = siftmatch(desc_1, desc_2);
matches_13 = siftmatch(desc_1, desc_3);
matches_14 = siftmatch(desc_1, desc_4);
p1 = [points_1(1:2, matches_12(1,:)); ones(1, length(matches_12))];
p2 = [points_2(1:2, matches_12(2,:)); ones(1, length(matches_12))];
[F_12, ~] = ransac_fundamental_matrix(p1, p2, 2.0, 1000); 
p1 = [points_1(1:2, matches_13(1,:)); ones(1, length(matches_13))];
p2 = [points_3(1:2, matches_13(2,:)); ones(1, length(matches_13))];
[F_13, ~] = ransac_fundamental_matrix(p1, p2, 2.0, 1000);
p1 = [points_1(1:2, matches_14(1,:)); ones(1, length(matches_14))];
p2 = [points_4(1:2, matches_14(2,:)); ones(1, length(matches_14))];
[F_14, ~] = ransac_fundamental_matrix(p1, p2, 2.0, 1000);


%% Plot the car trajectory (keypoint idx_car_I1 in image 1)

% ToDo: complete the code

idx_car_I1 = 1197;
idx_car_I2 = matches_12(2,matches_12(1,:)==idx_car_I1); % ToDo: identify the corresponding point of idx_car_I1 in image 2
idx_car_I3 = matches_13(2,matches_13(1,:)==idx_car_I1); % ToDo: identify the corresponding point of idx_car_I1 in image 3
idx_car_I4 = matches_14(2,matches_14(1,:)==idx_car_I1); % ToDo: identify the corresponding point of idx_car_I1 in image 4

% coordinates (in image 1) of the keypoint idx_car_I1 (point in a van). 
% point1_1 is the projection of a 3D point in the 3D trajectory of the van
point1_1 = [points_1(1:2,idx_car_I1)' 1]';
% coordinates (in image 1) of another 3D point in the same 3D trajectory of
% the van
point1_2 = [334 697 1]'; % (this is a given data)

% l1 is the projection of the 3D trajectory of keypoint idx_car_I1
% (it is the line that joins point1_1 and point1_2)
l1 = cross(point1_1, point1_2);% ToDo: compute the line
% plot the line
figure;imshow(im1);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(points_1(1,1197), points_1(2,1197), 'y*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 2
point2 = [points_2(1:2,idx_car_I2); 1];
% ToDo: compute the epipolar line of point2 in the reference image
l2 = F_12' * point2;
% plot the epipolar line
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'c');
% ToDo: compute the projection of point idx_car_I2 in the reference image 
pi2 = cross(l1,l2);
% plot this point
plot(pi2(1)/pi2(3), pi2(2)/pi2(3), 'c*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 3
point3 = [points_3(1:2,idx_car_I3); 1];
% ToDo: compute the epipolar line of point3 in the reference image
epipolar3 = F_13' * point3;
% plot the epipolar line
plot(t, -(epipolar3(1)*t + epipolar3(3)) / epipolar3(2), 'b');
% ToDo: compute the projection of point idx_car_I3 in the reference image
point3_im1 = cross(l1,epipolar3);
plot(point3_im1(1)/point3_im1(3), point3_im1(2)/point3_im1(3), 'b*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 4
point4 = [points_4(1:2,idx_car_I4); 1];
% ToDo: compute the epipolar line of point4 in the reference image
l4 = F_14' * point4;
% plot the epipolar line
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'g');
% ToDo: compute the projection of point idx_car_I4 in the reference image
point4_im1 = cross(l1,l4);
plot(point4_im1(1)/point4_im1(3), point4_im1(2)/point4_im1(3), 'g*');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. OPTIONAL: Photo-sequencing with your own images

% 4.1 Take a set of images of a moving scene from different viewpoints at 
%     different time instants. At least two images have to be taken from
%     roughly the same location by the same camera.

clear all;

% Read images & Reading unsorted images
im1rgb = imresize(imread('Data/IMG_0174.JPG'),[720 NaN]);
im2rgb = imresize(imread('Data/IMG_0175.JPG'),[720 NaN]);
im3rgb = imresize(imread('Data/IMG_0177.JPG'),[720 NaN]);
im4rgb = imresize(imread('Data/IMG_0176.JPG'),[720 NaN]);

im1 = sum(double(im1rgb), 3) / 3 / 255;
im2 = sum(double(im2rgb), 3) / 3 / 255;
im3 = sum(double(im3rgb), 3) / 3 / 255;
im4 = sum(double(im4rgb), 3) / 3 / 255;

% show images
figure;
subplot(2,2,1); imshow(im1rgb); axis image; title('Image 1');
subplot(2,2,2); imshow(im2rgb); axis image; title('Image 2');
subplot(2,2,3); imshow(im3rgb); axis image; title('Image 3');
subplot(2,2,4); imshow(im4rgb); axis image; title('Image 4');

% Compute SIFT keypoints
[points_1, desc_1] = sift(im1, 'Threshold', 0.015); % Do not change this threshold!
[points_2, desc_2] = sift(im2, 'Threshold', 0.015);
[points_3, desc_3] = sift(im3, 'Threshold', 0.015);
[points_4, desc_4] = sift(im4, 'Threshold', 0.015);

% Plot SIFT keypoints on im1
% figure; imshow(im1rgb); hold on; plot(points_1(1,:),points_1(2,:),'c*');

% 4.2 Implement the first part (until line 16) of the Algorithm 1 of the 
%     Photo-sequencing paper with a selection of the detected dynamic
%     features. You may reuse the code generated for the previous question.
%

%% Step 1: Match im1 and im2 (the similar images)
matches_12 = siftmatch(desc_1, desc_2);
figure; plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches_12, 'Stacking', 'v');

p1 = [points_1(1:2, matches_12(1,:)); ones(1, length(matches_12))];
p2 = [points_2(1:2, matches_12(2,:)); ones(1, length(matches_12))];

% [F, inliers] = ransac_fundamental_matrix(p1, p2, 2.0, 1000); 

% figure; plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches_12(:,inliers), 'Stacking', 'v');

%% Step 2: Classify the keypoints between dynamic and static
threshold = 60; % thresholdUSB = 30
% out_threshold = 400;
p1 = [points_1(1:2, matches_12(1,:)); ones(1, length(matches_12))];
p2 = [points_2(1:2, matches_12(2,:)); ones(1, length(matches_12))];

% p1 = [points_1(1:2, matches_12(1,inliers)); ones(1, length(matches_12(:,inliers)))];
% p2 = [points_2(1:2, matches_12(2,inliers)); ones(1, length(matches_12(:,inliers)))];
distances = sqrt(sum((p1-p2).^2,1));

static_12 = matches_12(:,distances<=threshold );%& distances<out_threshold);
dynamic_12 = matches_12(:,distances>threshold );%& distances<out_threshold);

figure; imshow(im1rgb);
hold on; plot(points_1(1,static_12(1,:)),points_1(2,static_12(1,:)),'r*');
plot(points_1(1,dynamic_12(1,:)),points_1(2,dynamic_12(1,:)),'b*');
plot(points_2(1,dynamic_12(2,:)),points_2(2,dynamic_12(2,:)),'g*');

%% Step 4: Match im1 with the other images
matches_13 = siftmatch(desc_1, desc_3);
matches_14 = siftmatch(desc_1, desc_4);

%% Step 5: Get the dynamic keypoints
dynamic_13 = matches_13(:,ismember(matches_13(1,:),dynamic_12(1,:)));
dynamic_14 = matches_14(:,ismember(matches_14(1,:),dynamic_12(1,:)));

static_13 = matches_13(:,ismember(matches_13(1,:),static_12(1,:)));
static_14 = matches_14(:,ismember(matches_14(1,:),static_12(1,:)));

%% Step 6: Compute fundamental matrices
p1 = [points_1(1:2, static_12(1,:)); ones(1, length(static_12(1,:)))];
p2 = [points_2(1:2, static_12(2,:)); ones(1, length(static_12(1,:)))];
[F_12, ~] = ransac_fundamental_matrix(p1, p2, 2.0, 1000); 
p1 = [points_1(1:2, static_13(1,:)); ones(1, length(static_13(1,:)))];
p2 = [points_3(1:2, static_13(2,:)); ones(1, length(static_13(1,:)))];
[F_13, ~] = ransac_fundamental_matrix(p1, p2, 2.0, 1000);
p1 = [points_1(1:2, static_14(1,:)); ones(1, length(static_14(1,:)))];
p2 = [points_4(1:2, static_14(2,:)); ones(1, length(static_14(1,:)))];
[F_14, ~] = ransac_fundamental_matrix(p1, p2, 2.0, 1000);

%% Step 7: 

p1 = [points_1(1:2, dynamic_12(1,:)); ones(1, length(dynamic_12(1,:)))];
p2 = [points_2(1:2, dynamic_12(2,:)); ones(1, length(dynamic_12(1,:)))];
directions = cross(p1,p2);

alpha3 = [];
alpha4 = [];
for i=1:size(dynamic_12,2)
%     f=figure; imshow(im1rgb); hold on;
    
    d_i = dynamic_12(1,i);
    point1 = [points_1(1:2,d_i); 1];
    point2 = [points_2(1:2,dynamic_12(2,i)); 1];
    ltrajectory = directions(:,i);
    
%     t=1:0.1:1000;
%     plot(t, -(ltrajectory(1)*t + ltrajectory(3)) / ltrajectory(2), 'y');
%     plot(point1(1), point1(2), 'y*');
%     plot(point2(1), point2(2), 'y*');

    % im3
    if max(dynamic_13(1,:)==d_i)
        point3 = [points_3(1:2,dynamic_13(2,dynamic_13(1,:)==d_i)); 1];
        
        epipolar3 = F_13' * point3;
        
        point3_im1 = cross(ltrajectory,epipolar3);
        point3_im1 = point3_im1 ./ point3_im1(3);
        
        alpha3 = [alpha3; sum((point3_im1-point1).^2)];
        
%         plot(t, -(epipolar3(1)*t + epipolar3(3)) / epipolar3(2), 'c');
%         plot(point3_im1(1), point3_im1(2), 'c*');
    end
    
    % im4
    if max(dynamic_14(1,:)==d_i)
        point4 = [points_4(1:2,dynamic_14(2,dynamic_14(1,:)==d_i)); 1];
        
        epipolar4 = F_14' * point4;
        
        point4_im1 = cross(ltrajectory,epipolar4);
        point4_im1 = point4_im1 ./ point4_im1(3);

        
        alpha4 = [alpha4; sum((point4_im1-point1).^2)];
        
%         plot(t, -(epipolar4(1)*t + epipolar4(3)) / epipolar4(2), 'b');
%         plot(point4_im1(1), point4_im1(2), 'b*');
    end
%     pause; close(f);
end

if median(alpha3)<median(alpha4)
    disp('Order: image1, image2, image3, image4');
else
    disp('Order: image1, image2, image4, image3');
end
%%
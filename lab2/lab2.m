%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 2: Image mosaics

addpath('sift');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Compute image correspondences

%% Open images

imargb = imread('Data/llanes/llanes_a.jpg');
imbrgb = imread('Data/llanes/llanes_b.jpg');
imcrgb = imread('Data/llanes/llanes_c.jpg');

% imargb = imread('Data/castle_int/0016_s.png');
% imbrgb = imread('Data/castle_int/0015_s.png');
% imcrgb = imread('Data/castle_int/0014_s.png');

% imargb = imread('Data/aerial/site13/frame00000.png');
% imbrgb = imread('Data/aerial/site13/frame00002.png');
% imcrgb = imread('Data/aerial/site13/frame00003.png');

ima = sum(double(imargb), 3) / 3 / 255;
imb = sum(double(imbrgb), 3) / 3 / 255;
imc = sum(double(imcrgb), 3) / 3 / 255;

% imargb = double(imread('Data/aerial/site22/frame_00001.tif'));
% imbrgb = double(imread('Data/aerial/site22/frame_00018.tif'));
% imcrgb = double(imread('Data/aerial/site22/frame_00030.tif'));
% ima = imargb;
% imb = imbrgb;
% imc = imcrgb;

%% Compute SIFT keypoints
[points_a, desc_a] = sift(ima, 'Threshold', 0.01);
[points_b, desc_b] = sift(imb, 'Threshold', 0.01);
[points_c, desc_c] = sift(imc, 'Threshold', 0.01);

figure;
imshow(imargb);%image(imargb)
hold on;
plot(points_a(1,:), points_a(2,:),'+y');
figure;
imshow(imbrgb);%image(imbrgb);
hold on;
plot(points_b(1,:), points_b(2,:),'+y');
figure;
imshow(imcrgb);%image(imcrgb);
hold on;
plot(points_c(1,:), points_c(2,:),'+y');

%% Match SIFT keypoints 

% between a and b
matches_ab = siftmatch(desc_a, desc_b);
figure;
plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), matches_ab, 'Stacking', 'v');

% between b and c
matches_bc = siftmatch(desc_b, desc_c);
figure;
plotmatches(imb, imc, points_b(1:2,:), points_c(1:2,:), matches_bc, 'Stacking', 'v');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Compute the homography (DLT algorithm) between image pairs

%% Compute homography (normalized DLT) between a and b, play with the homography
th = 3;
xab_a = [points_a(1:2, matches_ab(1,:)); ones(1, length(matches_ab))];
xab_b = [points_b(1:2, matches_ab(2,:)); ones(1, length(matches_ab))];
[Hab, inliers_ab] = ransac_homography_adaptive_loop(xab_a, xab_b, th, 1000); % ToDo: complete this function

figure;
plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), ...
    matches_ab(:,inliers_ab), 'Stacking', 'v');

vgg_gui_H(imargb, imbrgb, Hab);


%% Compute homography (normalized DLT) between b and c, play with the homography
xbc_b = [points_b(1:2, matches_bc(1,:)); ones(1, length(matches_bc))];
xbc_c = [points_c(1:2, matches_bc(2,:)); ones(1, length(matches_bc))];
[Hbc, inliers_bc] = ransac_homography_adaptive_loop(xbc_b, xbc_c, th, 1000); 

figure;
plotmatches(imb, imc, points_b(1:2,:), points_c(1:2,:), ...
    matches_bc(:,inliers_bc), 'Stacking', 'v');

vgg_gui_H(imbrgb, imcrgb, Hbc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Build the mosaic

corners = [-400 1200 -100 650];
iwb = apply_H_v2(imbrgb, eye(3), corners);   % ToDo: complete the call to the function
iwa = apply_H_v2(imargb, Hab, corners);    % ToDo: complete the call to the function
iwc = apply_H_v2(imcrgb, inv(Hbc), corners);    % ToDo: complete the call to the function

figure;
imshow(max(iwc, max(iwb, iwa)));%image(max(iwc, max(iwb, iwa)));axis off;
title('Mosaic A-B-C');

% ToDo: compute the mosaic with castle_int images
% ToDo: compute the mosaic with aerial images set 13
% ToDo: compute the mosaic with aerial images set 22
% ToDo: comment the results in every of the four cases: say why it works or
%       does not work

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Refine the homography with the Gold Standard algorithm

% Homography ab

x =  xab_a(1:2, inliers_ab); % ToDo: set the non-homogeneous point coordinates of the 
xp = xab_b(1:2, inliers_ab); % point correspondences we will refine with the geometric method
Xobs = [ x(:) ; xp(:) ];     % The column vector of observed values (x and x')
P0 = [ Hab(:) ; x(:) ];      % The parameters or independent variables

Y_initial = gs_errfunction( P0, Xobs ); % ToDo: create this function that we need to pass to the lsqnonlin function
% NOTE: gs_errfunction should return E(X) and not the sum-of-squares E=sum(E(X).^2)) that we want to minimize. 
% (E(X) is summed and squared implicitly in the lsqnonlin algorithm.) 
err_initial = sum( sum( Y_initial.^2 ));

options = optimset('Algorithm', 'levenberg-marquardt');
P = lsqnonlin(@(t) gs_errfunction(t, Xobs), P0, [], [], options);

Hab_r = reshape( P(1:9), 3, 3 );
f = gs_errfunction( P, Xobs ); % lsqnonlin does not return f
err_final = sum( sum( f.^2 ));

% we show the geometric error before and after the refinement
fprintf(1, 'Gold standard reproj error initial %f, final %f\n', err_initial, err_final);


%% See differences in the keypoint locations

% ToDo: compute the points xhat and xhatp which are the correspondences
% returned by the refinement with the Gold Standard algorithm

xhat = P(10:end);
length_xhat = size(xhat,1) / 2;
xhat = reshape(xhat, [2, length_xhat]);
xhat = [xhat; ones(1, length_xhat)];

xhatp = Hab_r * xhat;

figure;
imshow(imargb);image(imargb);
hold on;
plot(x(1,:), x(2,:),'+y');
plot(xhat(1,:), xhat(2,:),'+c');

figure;
imshow(imbrgb);image(imbrgb);
hold on;
plot(xp(1,:), xp(2,:),'+y');
plot(xhatp(1,:), xhatp(2,:),'+c');

%%  Homography bc

% ToDo: refine the homography bc with the Gold Standard algorithm
x =  xbc_b(1:2, inliers_bc); % ToDo: set the non-homogeneous point coordinates of the 
xp = xbc_c(1:2, inliers_bc); % point correspondences we will refine with the geometric method
Xobs = [ x(:) ; xp(:) ];     % The column vector of observed values (x and x')
P0 = [ Hbc(:) ; x(:) ];      % The parameters or independent variables

Y_initial = gs_errfunction( P0, Xobs ); % ToDo: create this function that we need to pass to the lsqnonlin function
% NOTE: gs_errfunction should return E(X) and not the sum-of-squares E=sum(E(X).^2)) that we want to minimize. 
% (E(X) is summed and squared implicitly in the lsqnonlin algorithm.) 
err_initial = sum( sum( Y_initial.^2 ));

options = optimset('Algorithm', 'levenberg-marquardt');
P = lsqnonlin(@(t) gs_errfunction(t, Xobs), P0, [], [], options);

Hbc_r = reshape( P(1:9), 3, 3 );
f = gs_errfunction( P, Xobs ); % lsqnonlin does not return f
err_final = sum( sum( f.^2 ));

% we show the geometric error before and after the refinement
fprintf(1, 'Gold standard reproj error initial %f, final %f\n', err_initial, err_final);

%% See differences in the keypoint locations

% ToDo: compute the points xhat and xhatp which are the correspondences
% returned by the refinement with the Gold Standard algorithm
xhat = P(10:end);
length_xhat = size(xhat,1) / 2;
xhat = reshape(xhat, [2, length_xhat]);
xhat = [xhat; ones(1, length_xhat)];

xhatp = Hbc_r * xhat;

figure;
imshow(imbrgb);%image(imbrgb);
hold on;
plot(x(1,:), x(2,:),'+y');
plot(xhat(1,:), xhat(2,:),'+c');

figure;
imshow(imcrgb);%image(imcrgb);
hold on;
plot(xp(1,:), xp(2,:),'+y');
plot(xhatp(1,:), xhatp(2,:),'+c');

%% Build mosaic
corners = [-400 1200 -100 650];
iwb = apply_H_v2(imbrgb, eye(3), corners); % ToDo: complete the call to the function
iwa = apply_H_v2(imargb, Hab_r, corners); % ToDo: complete the call to the function
iwc = apply_H_v2(imcrgb, inv(Hbc_r), corners); % ToDo: complete the call to the function

figure;
imshow(max(iwc, max(iwb, iwa)));%image(max(iwc, max(iwb, iwa)));axis off;
title('Mosaic A-B-C');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Calibration with a planar pattern

clear all;

%% Read template and images.
T     = imread('Data/calib/template.jpg');
I{1}  = imread('Data/calib/graffiti1.tif');
I{2}  = imread('Data/calib/graffiti2.tif');
I{3}  = imread('Data/calib/graffiti3.tif');
% I{4}  = imread('Data/calib/graffiti4.tif');
% I{5}  = imread('Data/calib/graffiti5.tif');
Tg = sum(double(T), 3) / 3 / 255;
Ig{1} = sum(double(I{1}), 3) / 3 / 255;
Ig{2} = sum(double(I{2}), 3) / 3 / 255;
Ig{3} = sum(double(I{3}), 3) / 3 / 255;

N = length(I);

%% Compute keypoints.
fprintf('Computing sift points in template... ');
[pointsT, descrT] = sift(Tg, 'Threshold', 0.05);
fprintf(' done\n');

points = cell(N,1);
descr = cell(N,1);
for i = 1:N
    fprintf('Computing sift points in image %d... ', i);
    [points{i}, descr{i}] = sift(Ig{i}, 'Threshold', 0.05);
    fprintf(' done\n');
end

%% Match and compute homographies.
H = cell(N,1);
for i = 1:N
    % Match against template descriptors.
    fprintf('Matching image %d... ', i);
    matches = siftmatch(descrT, descr{i});
    fprintf('done\n');

    % Fit homography and remove outliers.
    x1 = [pointsT(1:2, matches(1, :)); ones(1, length(matches))];
    x2 = [points{i}(1:2, matches(2, :)); ones(1, length(matches))];

    H{i} = 0;
    [H{i}, inliers] =  ransac_homography_adaptive_loop(x1, x2, 3, 1000);
%     H{i} = inv(H{i});
    % Plot inliers.
    figure;
    plotmatches(Tg, Ig{i}, pointsT(1:2,:), points{i}(1:2,:), matches(:, inliers));

    % Play with the homography
%     vgg_gui_H(T, I{i}, H{i});
end

%% Compute the Image of the Absolute Conic
v = [];
for i = 1:N
    h1 = H{i}(:,1);
    h2 = H{i}(:,2);
    h3 = H{i}(:,3);
    
    % vTij = (h1ih1j ; h1ih2j + h2ih1j ; h1ih3j + h3ih1j ; h2ih2j ; h2ih3j + h3ih2j ; h3ih3j)
    v12 = [h1(1)*h2(1), h1(1)*h2(2) + h1(2)*h2(1), h1(1)*h2(3)+h1(3)*h2(1),...
           h1(2)*h2(2), h1(2)*h2(3) + h1(3)*h2(2), h1(3)*h2(3)];
       
    v11 = [h1(1)*h1(1), h1(1)*h1(2) + h1(2)*h1(1), h1(1)*h1(3)+h1(3)*h1(1),...
           h1(2)*h1(2), h1(2)*h1(3) + h1(3)*h1(2), h1(3)*h1(3)];
       
    v22 = [h2(1)*h2(1), h2(1)*h2(2) + h2(2)*h2(1), h2(1)*h2(3)+h2(3)*h2(1),...
           h2(2)*h2(2), h2(2)*h2(3) + h2(3)*h2(2), h2(3)*h2(3)];
    v = [v; v12; v11 - v22];
end

[U, D, V] = svd(v); 
w = V(:,end);   

W = [w(1), w(2), w(3);
     w(2), w(4), w(5);
     w(3), w(5), w(6)]


%% Recover the camera calibration.

% Compute the matrix using the Cholesky factorization:
K = chol(inv(W),'upper');  % ToDo
    
% ToDo: in the report make some comments related to the obtained internal
%       camera parameters and also comment their relation to the image size

%% Compute camera position and orientation.
R = cell(N,1);
t = cell(N,1);
P = cell(N,1);
figure;hold;
for i = 1:N    
    h1 = H{i}(:,1);
    h2 = H{i}(:,2);
    h3 = H{i}(:,3);
    
    % ToDo: compute r1, r2, and t{i}
    r1 = K\h1;
    r2 = K\h2;
    t{i} = K\h3;
    
    % Solve the scale ambiguity by forcing r1 and r2 to be unit vectors.
    s = sqrt(norm(r1) * norm(r2)) * sign(t{i}(3));
    r1 = r1 / s;
    r2 = r2 / s;
    t{i} = t{i} / s;
    R{i} = [r1, r2, cross(r1,r2)];
    
    % Ensure R is a rotation matrix
    [U, S, V] = svd(R{i});
    R{i} = U * eye(3) * V';
   
    P{i} = K * [R{i} t{i}];
    plot_camera(P{i}, 800, 600, 200);
end

% ToDo: in the report explain how the optical center is computed in the
%       provided code

[ny,nx] = size(T);
p1 = [0 0 0]';
p2 = [nx 0 0]';
p3 = [nx ny 0]';
p4 = [0 ny 0]';
% Draw planar pattern
vgg_scatter_plot([p1 p2 p3 p4 p1], 'g');
% Paint image texture
surface('XData',[0 nx; 0 nx],'YData',[0 0; 0 0],'ZData',[0 0; -ny -ny],'CData',T,'FaceColor','texturemap');
colormap(gray);
axis equal;

%% Plot a static camera with moving calibration pattern.
figure; hold;
plot_camera(K * eye(3,4), 800, 600, 200);
% ToDo: complete the call to the following function with the proper
%       coordinates of the image corners in the new reference system

pts = [p1 p2 p3 p4 p1];
pts(end+1, : ) = 1;
for i = 1:N
%     Ti = eye(4);
%     Ti(1:3,end) = t{i};
%     Ri = [R{i}; [0 0 0]];
%     Ri = [Ri,[0 0 0 1]'];
%     ptsp = Ti*inv(Ri)*pts;
    ptsp = K\P{i}*pts
    vgg_scatter_plot( ptsp(1:3,:), 'r');
end

%% Augmented reality: Plot some 3D points on every camera.
[Th, Tw] = size(Tg);
cube = [0 0 0; 1 0 0; 1 0 0; 1 1 0; 1 1 0; 0 1 0; 0 1 0; 0 0 0; 0 0 1; 1 0 1; 1 0 1; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 0 1; 0 0 0; 1 0 0; 1 0 0; 1 0 1; 1 0 1; 0 0 1; 0 0 1; 0 0 0; 0 1 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 1 0; 0 0 0; 0 1 0; 0 1 0; 0 1 1; 0 1 1; 0 0 1; 0 0 1; 0 0 0; 1 0 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 1 0 1; 1 0 1; 1 0 0 ]';

X = (cube - .5) * Tw / 4 + repmat([Tw / 2; Th / 2; -Tw / 8], 1, length(cube));
X(end+1,:) = 1;

for i = 1:N
    figure; colormap(gray);
    imagesc(Ig{i});
    hold on;
    x = euclid(P{i} * X);
    vgg_scatter_plot(x, 'g');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. OPTIONAL: Detect the UPF logo in the two UPF images using the 
%%              DLT algorithm (folder "logos").
%%              Interpret and comment the results.

% Read images
imargb = imread('Data/logos/logoUPF.png');
imbrgb = imread('Data/logos/UPFbuilding.jpg');
imcrgb = imread('Data/logos/UPFstand.jpg');

ima = sum(double(imargb), 3) / 3 / 255;
imb = sum(double(imbrgb), 3) / 3 / 255;
imc = sum(double(imcrgb), 3) / 3 / 255;

% Compute SIFT keypoints
[points_a, desc_a] = sift(ima, 'Threshold', 0.01);
[points_b, desc_b] = sift(imb, 'Threshold', 0.01);
[points_c, desc_c] = sift(imc, 'Threshold', 0.01);

% Match SIFT keypoints 
matches_ab = siftmatch(desc_a, desc_b);
matches_ac = siftmatch(desc_a, desc_c);

% Compute homography
th = 3;
xab_a = [points_a(1:2, matches_ab(1,:)); ones(1, length(matches_ab))];
xab_b = [points_b(1:2, matches_ab(2,:)); ones(1, length(matches_ab))];
[Hab, inliers_ab] = ransac_homography_adaptive_loop(xab_a, xab_b, th, 1000);

xac_a = [points_a(1:2, matches_ac(1,:)); ones(1, length(matches_ac))];
xac_c = [points_c(1:2, matches_ac(2,:)); ones(1, length(matches_ac))];
[Hac, inliers_ac] = ransac_homography_adaptive_loop(xac_a, xac_c, th, 1000); 

% Compute logo corners
[ny,nx] = size(ima);
c1 = [0 0]';
c2 = [nx 0]';
c3 = [nx ny]';
c4 = [0 ny]';

corners = [c1 c2 c3 c4];
corners(end+1,:)=1;

cornersb = Hab*corners;
cornersb = cornersb./repmat(cornersb(end,:),3,1);
cornersc = Hac*corners;
cornersc = cornersc./repmat(cornersc(end,:),3,1);

% Plot logo corners
figure, imshow(imbrgb)
hold on,
plot(cornersb(1,:),cornersb(2,:),'k*');

figure, imshow(imcrgb)
hold on,
plot(cornersc(1,:),cornersc(2,:),'k*');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. OPTIONAL: Replace the logo of the UPF by the master logo
%%              in one of the previous images using the DLT algorithm.

logo = imread('Data/logos/logo_master.png');
[ny,nx,~] = size(logo);
c1 = [0 0]';
c2 = [nx 0]';
c3 = [nx ny]';
c4 = [0 ny]';

corners_master = [c1 c2 c3 c4];
corners_master(end+1,:)=1;

% Image b
Hb = homography2d(corners_master,cornersb);
corners = [1 size(imb,2) 1 size(imb,1)];
master_b = apply_H_v2(logo, Hb, corners);

figure;
imshow(max(master_b,imbrgb));

% Image c
Hc = homography2d(corners_master,cornersc);
corners = [1 size(imc,2) 1 size(imc,1)];
master_c = apply_H_v2(logo, Hc, corners);

figure;
imshow(max(master_c,imcrgb));

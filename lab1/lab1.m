%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 1: Image rectification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Applying image transformations

% ToDo: create the function  "apply_H" that gets as input a homography and
% an image and returns the image transformed by the homography.

% The size of the transformed image has to be automatically set so as to 
% contain the whole transformed image.

% At some point you will need to interpolate the image values at some points,
% you may use the Matlab function "interp2" for that.


%% 1.1. Similarities
I=imread('Data/0005_s.png'); % we have to be in the proper folder
% ToDo: generate a matrix H which produces a similarity transformation
% Parameters (4 dof)
theta = 0.1;    % rotation angle
s = 1;          % scaling factor
t = [1 1]';     % translation vector

R = [cos(theta) -sin(theta);
     sin(theta)  cos(theta)];
 
H = [s.*R, t;
     0, 0, 1];
 
I2 = apply_H(I, H);
figure(1);
subplot(1,2,1); imshow(I);         title('Original Image');
subplot(1,2,2); imshow(uint8(I2)); title('Similarity transformation');

%% 1.2. Affinities
% ToDo: generate a matrix H which produces an affine transformation
% Parameters (6 dof)
theta = 0.1;   %rotation angle
phi = 0.9;     %rotation angle
t = [1 1]';    % translation vector
d1 = 0.6;      % scaling factor
d2 = 0.9;      % scaling factor

R_theta = [cos(theta) -sin(theta);
           sin(theta)  cos(theta)]; 
R_phi = [cos(phi) -sin(phi);
         sin(phi)  cos(phi)]; 
D = [d1 0;
     0 d2];
 
A = R_theta * inv(R_phi) * D * R_phi;

H = [A, t;
     0, 0, 1];

I2 = apply_H(I, H);
%figure; imshow(I); figure; imshow(uint8(I2));
figure(2);
subplot(1,2,1); imshow(I);         title('Original Image');
subplot(1,2,2); imshow(uint8(I2)); title('Affine transformation');

%%
% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation
R_theta = [cos(theta) -sin(theta) 0;    % Rotation matrix 1
           sin(theta)  cos(theta) 0;
               0            0     1]; 
R_phi = [cos(phi) -sin(phi) 0;          % Rotation matrix 2
         sin(phi)  cos(phi) 0;
            0         0     1]; 
S = [d1 0 0;    % Scale matrix
     0 d2 0;
     0  0 1]; 
T = [1 0 t(1);  % Translation matrix
     0 1 t(2);
     0 0  1  ];
 
H_p = T * R_theta * inv(R_phi) * S * R_phi;

% ToDo: verify that the product of the four previous transformations
% produces the same matrix H as above

if(H == H_p)
    disp('H and H_p are equal')
else
    disp('H and H_p are diferent')
end

% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before

I2_p = apply_H(I, H_p);

figure;
subplot(1,2,1); imshow(uint8(I2)); title('I2');
subplot(1,2,2); imshow(uint8(I2_p)); title('I2\_p');

if(I2 == I2_p)
    disp('I2 and I2_p are equal')
else
    disp('I2 and I2_p are diferent')
end

%% 1.3 Projective transformations (homographies)
% ToDo: generate a matrix H which produces a projective transformation
% Parameters (8 dof)
theta = 0.1;   %rotation angle
phi = 0.1;     %rotation angle
t = [1 1]';    % translation vector
d1 = 1.5;      % scaling factor
d2 = 1;        % scaling factor
v = [0.0028 0.00028]';    % v

R_theta = [cos(theta) -sin(theta);
           sin(theta)  cos(theta)]; 
R_phi = [cos(phi) -sin(phi);
         sin(phi)  cos(phi)]; 
D = [d1 0;
     0 d2]; 
A = R_theta * inv(R_phi) * D * R_phi;

H = [A, t;
     v', 1];
 
I2 = apply_H(I, H);
%figure; imshow(I); figure; imshow(uint8(I2));
figure(1);
subplot(1,2,1); imshow(I);         title('Original Image');
subplot(1,2,2); imshow(uint8(I2)); title('Projective transformation');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification

% choose the image points
I=imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% indices of lines
i = 424;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 240;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 712;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 565;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = cross(p1,p2);
l1 = l1./l1(3);
l2 = cross(p3,p4);
l2 = l2./l2(3);
l3 = cross(p5,p6);
l3 = l3./l3(3);
l4 = cross(p7,p8);
l4 = l4./l4(3);

% show the chosen lines in the image
% figure;imshow(I);
% hold on;
% t=1:0.1:1000;
% plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
% plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
% plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
% plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

% ToDo: compute the homography that affinely rectifies the image
v1 = cross(l1,l2);
v2 = cross(l3,l4);

l_inf = cross(v1,v2);
l_inf = l_inf./l_inf(3);

H = [1 0 0; 0 1 0; l_inf(1) l_inf(2) l_inf(3)];

I2 = apply_H(I, H);%inv(H));
% figure; imshow(uint8(I2));

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
lr1 = inv(H)'*l1;
lr1 = lr1./lr1(3);

lr2 = inv(H)'*l2;
lr2 = lr2./lr2(3);

lr3 = inv(H)'*l3;
lr3 = lr3./lr3(3);

lr4 = inv(H)'*l4;
lr4 = lr4./lr4(3);

% show the transformed lines in the transformed image
figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation

% omega_inf = eye(3); omega_inf(3,3)=0;
% angle12 = acos(dot(omega_inf*l1,l2) / ...
%     ( sqrt(dot(omega_inf*l1,l1)) * sqrt(dot(omega_inf*l2,l2)) ) )*(180/pi);
% angle34 = acos(dot(omega_inf*l3,l4) / ...
%     ( sqrt(dot(omega_inf*l3,l3)) * sqrt(dot(omega_inf*l4,l4)) ) )*(180/pi);

angle12 = compute_angle(l1,l2);
angle34 = compute_angle(l3,l4);
fprintf('Angles before affine rectification: %.4f and %.4f\n',angle12,angle34);

% angle12p = acos(dot(omega_inf*lr1,lr2) / ...
%     ( sqrt(dot(omega_inf*lr1,lr1)) * sqrt(dot(omega_inf*lr2,lr2)) ) )*(180/pi);
% angle34p = acos(dot(omega_inf*lr3,lr4) / ...
%     ( sqrt(dot(omega_inf*lr3,lr3)) * sqrt(dot(omega_inf*lr4,lr4)) ) )*(180/pi);

angle12p = compute_angle(lr1,lr2);
angle34p = compute_angle(lr3,lr4);
fprintf('Angles after affine rectification: %.4f and %.4f\n',angle12p,angle34p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rectification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.

% Get the points in the corners of the window:
x1 = cross(lr1, lr3);
x2 = cross(lr2, lr4);
x3 = cross(lr1, lr4);
x4 = cross(lr2, lr3);

% With these points, compute the diagonal lines:
d1 = cross(x1, x2);
d2 = cross(x3, x4);

% Show the image and lines before metric rectification:
% figure;imshow(uint8(I2));
% hold on;
% t=1:0.1:1000;
% plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'g');
% plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'g');
% plot(t, -(d1(1)*t + d1(3)) / d1(2), 'y');
% plot(t, -(d2(1)*t + d2(3)) / d2(2), 'y');

% Matrix for the system of two equations:
B = [lr1(1) * lr3(1), lr1(1) * lr3(2) + lr1(2) * lr3(1), lr1(2) * lr3(2);
     d1(1) * d2(1), d1(1) * d2(2) + d1(2) * d2(1), d1(2) * d2(2)];

% Null vector/space of B obtained from the SVD.
s = null(B); 
S = [s(1), s(2);
     s(2), s(3)];
 
% Compute the upper triangular matrix using the Cholesky factorization:
K = chol(S, 'upper');

% Matrix for the metric rectification:
H_a = zeros(3,3);
H_a(1:2, 1:2) = K';
H_a(3,3) = 1;
H_a = inv(H_a); %Important

% Show the transformed image and lines:
I3 = apply_H(I2, H_a);

% Rectified lines, now with only a similiratiy distortion with the real
% world ones:
lm1 = inv(H_a)' * lr1;
lm1 = lm1./lm1(3);
lm2 = inv(H_a)' * lr2;
lm2 = lm2./lm2(3);
lm3 = inv(H_a)' * lr3;
lm3 = lm3./lm3(3);
lm4 = inv(H_a)' * lr4;
lm4 = lm4./lm4(3);
dm1 = inv(H_a)' * d1;
dm1 = dm1./dm1(3);
dm2 = inv(H_a)' * d2;
dm2 = dm2./dm2(3);

figure;imshow(uint8(I3));
hold on;
t=1:0.1:1000;
plot(t, -(lm1(1)*t + lm1(3)) / lm1(2), 'g');
plot(t, -(lm3(1)*t + lm3(3)) / lm3(2), 'g');
plot(t, -(lm2(1)*t + lm2(3)) / lm2(2), 'g');
plot(t, -(lm4(1)*t + lm4(3)) / lm4(2), 'g');
plot(t, -(dm1(1)*t + dm1(3)) / dm1(2), 'y');
plot(t, -(dm2(1)*t + dm2(3)) / dm2(2), 'y');

% Angles of lines before transformation:
ang_l1_l3_before = compute_angle(lr1, lr3);
ang_l2_l4_before = compute_angle(lr2, lr4);
ang_d1_d2_before = compute_angle(d1, d2);
fprintf('Angle of l1 and l3 before metric rectification = %f \n', ang_l1_l3_before)
fprintf('Angle of l2 and l4 before metric rectification = %f \n', ang_l2_l4_before)
fprintf('Angle of d1 and d2 before metric rectification = %f \n', ang_d1_d2_before)

% Angles of lines after transformation:
ang_l1_l3_after = compute_angle(lm1, lm3);
ang_l2_l4_after = compute_angle(lm2, lm4);
ang_d1_d2_after = compute_angle(dm1, dm2);
fprintf('Angle of l1 and l3 after metric rectification = %f \n', ang_l1_l3_after)
fprintf('Angle of l2 and l4 after metric rectification = %f \n', ang_l2_l4_after)
fprintf('Angle of d1 and d2 after metric rectification = %f \n', ang_d1_d2_after)

% [I3, H_a] = metricRectification(I2, H, A, 424, 240, 712, 565);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Affine and Metric Rectification of the left facade of image 0001

% ToDo: Write the code that rectifies the left facade of image 0001 with
%       the stratified method (affine + metric). 
%       Crop the initial image so that only the left facade is visible.
%       Show the (properly) transformed lines that use in every step.

% choose the image points
I = imread('Data/0001_s.png');
A = load('Data/0001_s_info_lines.txt');
% indices of lines
i = 614;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 159;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 645;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 541;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

% AFFINE
% Compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = cross(p1,p2);
l1 = l1./l1(3);
l2 = cross(p3,p4);
l2 = l2./l2(3);
l3 = cross(p5,p6);
l3 = l3./l3(3);
l4 = cross(p7,p8);
l4 = l4./l4(3);

% show the chosen lines in the image
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

% Compute the homography that affinely rectifies the image
v1 = cross(l1,l2);
v2 = cross(l3,l4);

l_inf = cross(v1,v2);
l_inf = l_inf./l_inf(3);

H = [1 0 0; 0 1 0; l_inf(1) l_inf(2) l_inf(3)];

I2 = apply_H(I, H);
% figure; imshow(uint8(I2));

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
lr1 = inv(H)'*l1;
lr1 = lr1./lr1(3);

lr2 = inv(H)'*l2;
lr2 = lr2./lr2(3);

lr3 = inv(H)'*l3;
lr3 = lr3./lr3(3);

lr4 = inv(H)'*l4;
lr4 = lr4./lr4(3);

% show the transformed lines in the transformed image
figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');

% Compute angles
angle12 = compute_angle(l1,l2);
angle34 = compute_angle(l3,l4);
angle12p = compute_angle(lr1,lr2);
angle34p = compute_angle(lr3,lr4);
fprintf('Angles before affine rectification: %.4f and %.4f\n',angle12,angle34);
fprintf('Angles after affine rectification: %.4f and %.4f\n',angle12p,angle34p);

%METRIC
% Get the points in the corners of the window:
x1 = cross(lr1, lr3);
x2 = cross(lr2, lr4);
x3 = cross(lr1, lr4);
x4 = cross(lr2, lr3);

% With these points, compute the diagonal lines:
d1 = cross(x1, x2);
d2 = cross(x3, x4);

% Show the image and lines before metric rectification:
figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'g');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'g');
plot(t, -(d1(1)*t + d1(3)) / d1(2), 'y');
plot(t, -(d2(1)*t + d2(3)) / d2(2), 'y');

% Matrix for the system of two equations:
B = [lr1(1) * lr3(1), lr1(1) * lr3(2) + lr1(2) * lr3(1), lr1(2) * lr3(2);
     d1(1) * d2(1), d1(1) * d2(2) + d1(2) * d2(1), d1(2) * d2(2)];

% Null vector/space of B obtained from the SVD.
s = null(B); 
S = [s(1), s(2);
     s(2), s(3)];
 
% Compute the upper triangular matrix using the Cholesky factorization:
K = chol(S, 'upper');

% Matrix for the metric rectification:
H_a = zeros(3,3);
H_a(1:2, 1:2) = K';
H_a(3,3) = 1;
H_a = inv(H_a); %Important

% Show the transformed image and lines:
I3 = apply_H_2(I2, H_a,transformed_corners2);

% Rectified lines, now with only a similiratiy distortion with the real
% world ones:
lm1 = inv(H_a)' * lr1;
lm1 = lm1./lm1(3);
lm2 = inv(H_a)' * lr2;
lm2 = lm2./lm2(3);
lm3 = inv(H_a)' * lr3;
lm3 = lm3./lm3(3);
lm4 = inv(H_a)' * lr4;
lm4 = lm4./lm4(3);
dm1 = inv(H_a)' * d1;
dm1 = dm1./dm1(3);
dm2 = inv(H_a)' * d2;
dm2 = dm2./dm2(3);

figure;imshow(uint8(I3));
hold on;
t=1:0.1:1000;
plot(t, -(lm1(1)*t + lm1(3)) / lm1(2), 'g');
plot(t, -(lm3(1)*t + lm3(3)) / lm3(2), 'g');
plot(t, -(lm2(1)*t + lm2(3)) / lm2(2), 'g');
plot(t, -(lm4(1)*t + lm4(3)) / lm4(2), 'g');
plot(t, -(dm1(1)*t + dm1(3)) / dm1(2), 'y');
plot(t, -(dm2(1)*t + dm2(3)) / dm2(2), 'y');

% Angles of lines before transformation:
ang_l1_l3_before = compute_angle(lr1, lr3);
ang_l2_l4_before = compute_angle(lr2, lr4);
ang_d1_d2_before = compute_angle(d1, d2);
fprintf('Angle of l1 and l3 before metric rectification = %f \n', ang_l1_l3_before)
fprintf('Angle of l2 and l4 before metric rectification = %f \n', ang_l2_l4_before)
fprintf('Angle of d1 and d2 before metric rectification = %f \n', ang_d1_d2_before)

% Angles of lines after transformation:
ang_l1_l3_after = compute_angle(lm1, lm3);
ang_l2_l4_after = compute_angle(lm2, lm4);
ang_d1_d2_after = compute_angle(dm1, dm2);
fprintf('Angle of l1 and l3 after metric rectification = %f \n', ang_l1_l3_after)
fprintf('Angle of l2 and l4 after metric rectification = %f \n', ang_l2_l4_after)
fprintf('Angle of d1 and d2 after metric rectification = %f \n', ang_d1_d2_after)

%%
corner1 = H_a*[1 1 1]';
corner2 = H_a*[1 size(I2,1) 1]';

corner2(2)-corner1(2)

%%
[rows, cols, chan] = size(I2);
corners = zeros(4,3);
corners(1,:) = [1 1 1]';
corners(2,:) = [cols 1 1]';
corners(3,:) = [1 rows 1]';
corners(4,:) = [cols rows 1]';


% Transform the image corners 
transformed_corners = H*corners';

% Make transformed_corners homogeneous
for i=1:4,
    transformed_corners(:,i) = transformed_corners(:,i)/transformed_corners(3,i);
end

% Compute the new extrenal transformed corner coordinates      
xmin = round(min([transformed_corners(1,1), transformed_corners(1,2), ...
                  transformed_corners(1,3), transformed_corners(1,4)]));

xmax = round(max([transformed_corners(1,1), transformed_corners(1,2), ...
                  transformed_corners(1,3), transformed_corners(1,4)]));

ymin = round(min([transformed_corners(2,1), transformed_corners(2,2), ...
                  transformed_corners(2,3), transformed_corners(2,4)]));
ymax = round(max([transformed_corners(2,1), transformed_corners(2,2), ...
                  transformed_corners(2,3), transformed_corners(2,4)]));
                  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)

% choose the image points
I=imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% indices of lines
i = 424;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 240;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 712;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 565;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

i = 227;
p9 = [A(i,1) A(i,2) 1]';
p10 = [A(i,3) A(i,4) 1]';
i = 367;
p11 = [A(i,1) A(i,2) 1]';
p12 = [A(i,3) A(i,4) 1]';
i = 534;
p13 = [A(i,1) A(i,2) 1]';
p14 = [A(i,3) A(i,4) 1]';
i = 576;
p15 = [A(i,1) A(i,2) 1]';
p16 = [A(i,3) A(i,4) 1]';

% Compute the lines 
l1 = cross(p1,p2);
% l1 = l1./l1(3);
l2 = cross(p3,p4);
% l2 = l2./l2(3);
l3 = cross(p5,p6);
% l3 = l3./l3(3);
l4 = cross(p7,p8);
% l4 = l4./l4(3);

x1 = cross(l1,l3);
x1 = x1./x1(3);
x2 = cross(l2,l4);
x2 = x2./x2(3);
x3 = cross(l1,l4);
x3 = x3./x3(3);
x4 = cross(l2,l3);
x4 = x4./x4(3);

d1 = cross(x1, x2);
d2 = cross(x3, x4);

l5 = cross(p9,p10);
% l5 = l5./l5(3);
l6 = cross(p11,p12);
% l6 = l6./l6(3);
l7 = cross(p13,p14);
% l7 = l7./l7(3);
l8 = cross(p15,p16);
% l8 = l8./l8(3);

x5 = cross(l5,l7);
x5 = x5./x5(3);
x6 = cross(l6,l8);
x6 = x6./x6(3);
x7 = cross(l5,l8);
x7 = x7./x7(3);
x8 = cross(l6,l7);
x8 = x8./x8(3);

d3 = cross(x5, x6);
d4 = cross(x7, x8);

% show the chosen lines in the image
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');
plot(t, -(l5(1)*t + l5(3)) / l5(2), 'y');
plot(t, -(l6(1)*t + l6(3)) / l6(2), 'y');
plot(t, -(l7(1)*t + l7(3)) / l7(2), 'y');
plot(t, -(l8(1)*t + l8(3)) / l8(2), 'y');
plot(t, -(d1(1)*t + d1(3)) / d1(2), 'g');
plot(t, -(d2(1)*t + d2(3)) / d2(2), 'g');
plot(t, -(d3(1)*t + d3(3)) / d3(2), 'g');
plot(t, -(d4(1)*t + d4(3)) / d4(2), 'g');

B = [
    %l1 l3
    l1(1)*l3(1),...
    (l1(1)*l3(2) + l1(2)*l3(1))/2,...
    l1(2)*l3(2),...
    (l1(1)*l3(3)+ l1(3)*l3(1))/2,...
    (l1(2)*l3(3)+ l1(3)*l3(2))/2,...
    l1(3)*l3(3);
    
    %d1 d2
    d1(1)*d2(1),...
    (d1(1)*d2(2) + d1(2)*d2(1))/2,...
    d1(2)*d2(2),...
    (d1(1)*d2(3)+ d1(3)*d2(1))/2,...
    (d1(2)*d2(3)+ d1(3)*d2(2))/2,...
    d1(3)*d2(3);
    
    %l5 l7
    l5(1)*l7(1),...
    (l5(1)*l7(2) + l5(2)*l7(1))/2,...
    l5(2)*l7(2),...
    (l5(1)*l7(3)+ l5(3)*l7(1))/2,...
    (l5(2)*l7(3)+ l5(3)*l7(2))/2,...
    l5(3)*l7(3);
    
    %d3 d4
    d3(1)*d4(1),...
    (d3(1)*d4(2) + d3(2)*d4(1))/2,...
    d3(2)*d4(2),...
    (d3(1)*d4(3)+ d3(3)*d4(1))/2,...
    (d3(2)*d4(3)+ d3(3)*d4(2))/2,...
    d3(3)*d4(3);
    
    %l2 l4
    l2(1)*l4(1),...
    (l2(1)*l4(2) + l2(2)*l4(1))/2,...
    l2(2)*l4(2),...
    (l2(1)*l4(3)+ l2(3)*l4(1))/2,...
    (l2(2)*l4(3)+ l2(3)*l4(2))/2,...
    l2(3)*l4(3)];
 
%(l1m1,(l1m2 + l2m1)/2, l2m2,(l1m3 + l3m1)/2,(l2m3 + l3m2)/2, l3m3) c = 0

c = null(B); % Null vector/space of B obtained from the SVD.

% C =
% a b/2 d/2
% b/2 c e/2
% d/2 e/2 f

C = [c(1), c(2)/2, c(4)/2;
    c(2)/2, c(3), c(5)/2;
    c(4)/2, c(5), c(6)];

%%
[U,S,V] = svd(C);

H = V';
H = H./H(3,3);
I2 = apply_H(I,H);

figure,imshow(uint8(I2));



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
figure; imshow(I); figure; imshow(uint8(I2));

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
figure; imshow(I); figure; imshow(uint8(I2));
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
figure; imshow(I); figure; imshow(uint8(I2));

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
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

% ToDo: compute the homography that affinely rectifies the image
v1 = cross(l1,l2);
v2 = cross(l3,l4);

l_inf = cross(v1,v2);
% l_inf = l_inf./sqrt(l_inf(1)^2+l_inf(2)^2);
l_inf = l_inf./l_inf(3);

H = [1 0 0; 0 1 0; l_inf(1) l_inf(2) l_inf(3)];

I2 = apply_H(I, H);%inv(H));
figure; imshow(uint8(I2));

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
lr1 = inv(H)'*l1;
% lr1./sqrt(lr1(1)^2+lr1(2)^2)
lr1 = lr1./lr1(3);

lr2 = inv(H)'*l2;
% lr2./sqrt(lr2(1)^2+lr2(2)^2)
lr2 = lr2./lr2(3);

lr3 = inv(H)'*l3;
% lr3./sqrt(lr3(1)^2+lr3(2)^2)
lr3 = lr3./lr3(3);

lr4 = inv(H)'*l4;
% lr4./sqrt(lr4(1)^2+lr4(2)^2)
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

%Angles before transformation
angle12 = acos(l1(1)*l2(1) + l1(2)*l2(2))*(180/pi)
angle34 = acos(l3(1)*l4(1) + l3(2)*l4(2))*(180/pi)

%Angles after transformation
angle12p = acos(lr1(1)*lr2(1) + lr1(2)*lr2(2))*(180/pi)
angle34p = acos(lr3(1)*lr4(1) + lr3(2)*lr4(2))*(180/pi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rectification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Affine and Metric Rectification of the left facade of image 0001

% ToDo: Write the code that rectifies the left facade of image 0001 with
%       the stratified method (affine + metric). 
%       Crop the initial image so that only the left facade is visible.
%       Show the (properly) transformed lines that use in every step.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)




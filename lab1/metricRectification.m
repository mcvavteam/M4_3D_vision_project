% This function gets as input an affine rectified image (I2) with
% the projectivity that performs this rectification (H_p), and a set 
% of perpendicular points
function [I3, H_a] = metricRectification(I2, H_p, A, i1, i2, i3, i4)
% indices of lines
p1 = [A(i1,1) A(i1,2) 1]';
p2 = [A(i1,3) A(i1,4) 1]';
p3 = [A(i2,1) A(i2,2) 1]';
p4 = [A(i2,3) A(i2,4) 1]';
p5 = [A(i3,1) A(i3,2) 1]';
p6 = [A(i3,3) A(i3,4) 1]';
p7 = [A(i4,1) A(i4,2) 1]';
p8 = [A(i4,3) A(i4,4) 1]';

% compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = cross(p1,p2);
l2 = cross(p3,p4);
l3 = cross(p5,p6);
l4 = cross(p7,p8);

% Transform the lines (affine):
l1_a = inv(H_p)' * l1;
l1_a = l1_a./l1_a(3);
l2_a = inv(H_p)' * l2;
l2_a = l2_a./l2_a(3);
l3_a = inv(H_p)' * l3;
l3_a = l3_a./l3_a(3);
l4_a = inv(H_p)' * l4;
l4_a = l4_a./l4_a(3);

figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(l1_a(1)*t + l1_a(3)) / l1_a(2), 'y');
plot(t, -(l2_a(1)*t + l2_a(3)) / l2_a(2), 'y');
plot(t, -(l3_a(1)*t + l3_a(3)) / l3_a(2), 'y');
plot(t, -(l4_a(1)*t + l4_a(3)) / l4_a(2), 'y');

% Get the points in the corners of the window:
x1 = cross(l1_a, l3_a);
x2 = cross(l2_a, l4_a);
x3 = cross(l1_a, l4_a);
x4 = cross(l2_a, l3_a);

% With these points, compute the diagonal lines:
d1 = cross(x1, x2);
d2 = cross(x3, x4);

% Show the image and lines before metric rectification:
figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(l1_a(1)*t + l1_a(3)) / l1_a(2), 'g');
plot(t, -(l3_a(1)*t + l3_a(3)) / l3_a(2), 'g');
plot(t, -(d1(1)*t + d1(3)) / d1(2), 'y');
plot(t, -(d2(1)*t + d2(3)) / d2(2), 'y');

% Matrix for the system of two equations:
B = [l1_a(1) * l3_a(1), l1_a(1) * l3_a(2) + l1_a(2) * l3_a(1), l1_a(2) * l3_a(2);
     d1(1) * d2(1), d1(1) * d2(2) + d1(2) * d2(1), d1(2) * d2(2)];

s = null(B); % Null vector/space of B obtained from the SVD.
S = [s(1), s(2);
     s(2), s(3)];

% Compute the upper triangular matrix using the Cholesky factorization:
% K = inv(chol(inv(S))); % Correcto
% K = inv(chol(S)); % Incorrecto
K = chol(S, 'upper');

% Matrix for the metric rectification:
% K of H_a should be invers but Idk peta.
H_a = zeros(3,3);
H_a(1:2, 1:2) = K';
H_a(3,3) = 1;

% Importante:
H_a = inv(H_a);

% Show the transformed image and lines:
I3 = apply_H(I2, H_a);

% Rectified lines, now with only a similiratiy distortion with the real
% world ones:
l1_s = inv(H_a)' * l1_a;
%l1_s = l1_s./l1_s(3);
l2_s = inv(H_a)' * l2_a;
%l2_s = l2_s./l2_s(3);
l3_s = inv(H_a)' * l3_a;
%l3_s = l3_s./l3_s(3);
l4_s = inv(H_a)' * l4_a;
%l4_s = l4_s./l4_s(3);
d1_s = inv(H_a)' * d1;
%d1_s = d1_s./d1_s(3);
d2_s = inv(H_a)' * d2;
%d2_s = d2_s./d2_s(3);


figure;imshow(uint8(I3));
hold on;
t=1:0.1:1000;
plot(t, -(l1_s(1)*t + l1_s(3)) / l1_s(2), 'g');
plot(t, -(l3_s(1)*t + l3_s(3)) / l3_s(2), 'g');
plot(t, -(d1_s(1)*t + d1_s(3)) / d1_s(2), 'y');
plot(t, -(d2_s(1)*t + d2_s(3)) / d2_s(2), 'y');

% Angles of lines before transformation:
ang_l1_l3_before = compute_angle(l1_a, l3_a);
ang_l2_l4_before = compute_angle(l2_a, l4_a);
ang_d1_d2_before = compute_angle(d1, d2);
fprintf('Angle of l1 and l3 before metric rectification = %f \n', ang_l1_l3_before)
fprintf('Angle of l2 and l4 before metric rectification = %f \n', ang_l2_l4_before)
fprintf('Angle of d1 and d2 before metric rectification = %f \n', ang_d1_d2_before)

% Angles of lines after transformation:
ang_l1_l3_after = compute_angle(l1_s, l3_s);
ang_l2_l4_after = compute_angle(l2_s, l4_s);
ang_d1_d2_after = compute_angle(d1_s, d2_s);
fprintf('Angle of l1 and l3 after metric rectification = %f \n', ang_l1_l3_after)
fprintf('Angle of l2 and l4 after metric rectification = %f \n', ang_l2_l4_after)
fprintf('Angle of d1 and d2 after metric rectification = %f \n', ang_d1_d2_after)

return
    
end
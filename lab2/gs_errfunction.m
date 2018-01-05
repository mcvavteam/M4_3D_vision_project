function [ Error ] = gs_errfunction( p0, Xobs )
% Arguments:
% - Xobs = [ x(:) ; xp(:) ];     % The column vector of observed values (x and x')
% - P0 = [ Hab(:) ; x(:) ];      % The parameters or independent variables

% get the homography
H = reshape( p0(1:9), [3,3]);

points_length = size(Xobs, 1) / 2;
x  = Xobs(1:points_length);
xp = Xobs(points_length+1:end);
x  = reshape(x, [2, size(x,1)/2]);
xp = reshape(xp, [2, size(xp,1)/2]);

% Get x hat
% the first 9 elements correspond to the homography
x_hat = p0(10:end);
length_xhat = size(x_hat, 1) / 2;
x_hat = reshape(x_hat, [2, length_xhat]);
x_hat = [x_hat; ones(1, length_xhat)];


% Calculate Homo*x_hat and the error between estimations and observations
x_hat_p = H * x_hat;

error1 = abs(x - euclid(x_hat));
error2 = abs(xp - euclid(x_hat_p));
Error = [error1(:) error2(:)];

end
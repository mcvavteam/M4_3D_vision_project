function [ Xproj, Pproj ] = factorization_method( x, init )
% computes a preconstruction with the factorization method of Sturm and
% Triggs '1996
%
% Arguments: x: Number of cameras (3xNPoint matrix each camera).
% This function returns an estimate of:
%       Pproj: 3*Ncam x 4 matrix containing the camera matrices
%       Xproj: 4 x Npoints matrix of homogeneous coordinates of 3D points
% 
% As a convergence criterion you may compute the Euclidean
% distance (d) between data points and projected points in both images 
% and stop when (abs(d - d_old)/d) < 0.1 where d_old is the distance
% in the previous iteration.

NumCams = length(x);
NumPoints = size(x{1},2);

% Normalize the points in each image (Similarity transformation.)
for i=1:NumCams
    [x_norm{i} H{i}] = normalise2dpts(x{i});
end

% Initialize lambda (put a flag to initialize with ones or Sturm)
lambda = ones(NumCams, NumPoints);
if isequal(init, 'sturm')
    % Initialize all projective depths 
    for i = 1:NumCams
        F{i} = fundamental_matrix( x{i}, x{1} );
        [U, D, V] = svd(F{i});
        e{i} = V(:,3) / V(3,3);
    end

    for i = 1:NumCams
        for j = 1:NumPoints
            numerator = x{1}(:, j)' * F{i} * cross(e{i}, x{i}(:,j));
            denominator = norm(cross(e{i}, x{i}(:,j))).^2 * lambda(1, j);
            lambda(i,j) = numerator / denominator;
        end
    end
end

% Initialize d and d_old
d_old = 9e+10;
d = 0;

% Initialize the measurement that multiplies lambda to obtain the
% measurement matrix M

x_meas = [];
   
for i = 1:NumCams
    x_meas = [x_meas; x_norm{i}];
end

% Iterate until teh convergence criterion 
tol = 0.1;
while double(abs(d - d_old)/d) > tol
    % Lecture: Alternate rescaling the rows for the depth matrix (lambda) to have unit norm and the columns of depth matrix to have unit norm until it stops changing significantly (2 loops)
    % Reescale Lambda 
    for i = 1:2
        for j = 1:size(lambda,1)
            lambda(j,:) = lambda(j,:) / norm(lambda(j,:));
        end
        for j = 1:size(lambda,2)
            lambda(:,j) = lambda(:,j) / norm(lambda(:,j));
        end
    end
    
    aux = [];
    for i = 1:NumCams
        aux = [aux; lambda(i,:); lambda(i,:); lambda(i,:)];
    end
    lambda = aux;
    d_old = d;
    % Measuremente Matrix M
    M = lambda.*x_meas;
    
    [U, D, V] = svd(M);
    
    % Now PM = UD4 and XM = V4T
    % Camera projection matrices and homogenous 3D points:
    Pproj = U * D(:,1:4);
    Xproj = V(:,1:4)';
    
    % If sum_i sum_j of d(xji , P^i X_j)2 converges then stop; 
    % Otherwise let lambda_ij = (P^iX_j)_3.
    % Compute distance between original 2D points and projected ones:
    d = sum(sum((x_meas - Pproj * Xproj).^2));
    
    P = {};
    P{1} = Pproj(1:3,:);
    P{2} = Pproj(4:6,:);
    lambda = [];
    for i = 1:length(P);
        aux = P{i} * Xproj;
        lambda = [lambda; aux(3,:)];
    end
end


% Unnormalize the camera matrices 
Pproj = [];
Pproj2 = [];
for i = 1:NumCams
    %Pproj = [Pproj; (H{i}\P{i})];
    Pproj = [Pproj; (inv(H{i})*P{i})];
end


end


function [F, idx_inliers] = ransac_fundamental_matrix(x1, x2, th, max_it)

[Ncoords, Npoints] = size(x1);

% ransac
it = 0;
best_inliers = [];
% probability that at least one random sample set is free of outliers
p = 0.999; 
while it < max_it
    
    points = randomsample(Npoints, 8);
    F = fundamental_matrix(x1(:,points), x2(:,points));
    inliers = compute_inliers(F, x1, x2, th);
    
    % test if it is the best model so far
    if length(inliers) > length(best_inliers)
        best_inliers = inliers;
    end    
    
    % update estimate of max_it (the number of trials) to ensure we pick, 
    % with probability p, an initial data set with no outliers
    fracinliers =  length(inliers)/Npoints;
    pNoOutliers = 1 -  fracinliers^4;
    pNoOutliers = max(eps, pNoOutliers);  % avoid division by -Inf
    pNoOutliers = min(1-eps, pNoOutliers);% avoid division by 0
    max_it = log(1-p)/log(pNoOutliers);
    
    it = it + 1;
end

% compute F for the best inlier
F = fundamental_matrix(x1(:,best_inliers), x2(:,best_inliers));
idx_inliers = best_inliers;

% Inliers are obtained with a threshold on the first order approx.
% of the geometi error (Sampson distance)
function idx_inliers = compute_inliers(F, x1, x2, th)
    % Check that H is invertible
    num_p = size(x1, 2);
       
    p2Fp1 = zeros(1, num_p);
    
    for i = 1: num_p
        p2Fp1(i) = x2(:,i)' * F * x1(:,i);
    end
    
    Fp1 = F*x1;
    Fp2 = F'*x2;
   
    % Computing Sampson Distance (Error)
    d = p2Fp1.^2 ./ (Fp1(1,:).^2 + Fp1(2,:).^2 + Fp2(1,:).^2 + Fp2(2,:).^2); 
    idx_inliers = find(d < th.^2);


% function xn = normalise(x)    
%     xn = x ./ repmat(x(end,:), size(x,1), 1);

    
function item = randomsample(npts, n)
	a = [1:npts]; 
    item = zeros(1,n);    
    for i = 1:n
	  % Generate random value in the appropriate range 
	  r = ceil((npts-i+1).*rand);
	  item(i) = a(r);       % Select the rth element from the list
	  a(r)    = a(end-i+1); % Overwrite selected element
    end                       % ... and repeat
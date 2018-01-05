function [H, idx_inliers] = ransac_homography_adaptive_loop(x1, x2, th, max_it)

    [Ncoords, Npoints] = size(x1);

    % ransac
    it = 0;
    best_inliers = [];
    while it < max_it

        points = randomsample(Npoints, 4);
        H = homography2d(x1(:,points), x2(:,points)); % ToDo: you have to create this function
        inliers = compute_inliers(H, x1, x2, th);

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
        p=0.99;
        max_it = log(1-p)/log(pNoOutliers);

        it = it + 1;
    end

    % compute H from all the inliers
    H = homography2d(x1(:,best_inliers), x2(:,best_inliers));
    idx_inliers = best_inliers;
end
%% compute_inliers
function idx_inliers = compute_inliers(H, x1, x2, th)
    % Check that H is invertible
    if abs(log(cond(H))) > 15
        idx_inliers = [];
        return
    end
    
    % compute the symmetric geometric error
    H_x1 = normalise(H * x1);
    % Reason of inv(H) inf. loop.
    H_x2 = normalise(inv(H) * x2);
    x1 = normalise(x1);
    x2 = normalise(x2);
    
    % ToDo
    % x2 to H_x1 and H_x2 to x1
    d2 = ((x2(1,:)-H_x1(1,:)).^2 + (x2(2,:)-H_x1(2,:)).^2 + (x2(3,:)-H_x1(3,:)).^2) ...
        + ((H_x2(1,:)-x1(1,:)).^2 + (H_x2(2,:)-x1(2,:)).^2 + (H_x2(3,:)-x1(3,:)).^2);
    %d2 = % ToDo
    idx_inliers = find(d2 < th.^2);
end
%% normalise
function xn = normalise(x)    
    xn = x ./ repmat(x(end,:), size(x,1), 1);
end
%% randomsample
function item = randomsample(npts, n)
	a = [1:npts]; 
    item = zeros(1,n);    
    for i = 1:n
	  % Generate random value in the appropriate range 
	  r = ceil((npts-i+1).*rand);
	  item(i) = a(r);       % Select the rth element from the list
	  a(r)    = a(end-i+1); % Overwrite selected element
    end                       % ... and repeat
end
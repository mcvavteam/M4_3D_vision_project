function [ F ] = fundamental_matrix(pt1, pt2)
% Fundamental matrix (8-point algorithm)

    % Number of points
    [~,n1] = size(pt1);
    [~,n2] = size(pt2);
    n = min(n1,n2);

    % Attempt to normalise each set of points so that the origin 
    % is at centroid and mean distance from origin is sqrt(2).
    [pt1, T1] = normalise2dpts(pt1);
    [pt2, T2] = normalise2dpts(pt2);

    % Construct W
    u1 = pt1(1,:)';
    u2 = pt2(1,:)';
    v1 = pt1(2,:)';
    v2 = pt2(2,:)';
    
    W = [u1.*u2, v1.*u2, u2, u1.*v2, v1.*v2, v2, u1, v1, ones(n,1)];

    [~,~,V] = svd(W);
    f = V(:,end);
    F = reshape(f,3,3)';
    
    % Change fo rank (from rank 3 to rank 2)
    [U,D,V] = svd(F);
    D(3,3)=0;    
    F = U*D*V';
    
    % De-normalize
    F = T2'*F*T1;
end


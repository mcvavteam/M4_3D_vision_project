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

    % Construct A
%     W = zeros(n,9);
        
    W = [];

    for i=1:n
        u =pt1(1,i);
        v =pt1(2,i);
        
        up =pt2(2,i);
        vp =pt2(2,i);
        
%         w = [uu' vu' u' uv' vv' v' u v 1]
        w = [u*up v*up up u*vp v*vp vp u v 1];
        W = [W;w];
    end

    [U,D,V] = svd(W);
    f = V(:,end);

    F = reshape(f,3,3)';
    F = T2\F*T1;
    F = F / F(3,3);

end


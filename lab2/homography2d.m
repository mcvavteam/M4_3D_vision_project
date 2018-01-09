function [ H ] = homography2d(pt1, pt2)
% Direct Linear Transform (DLT) algorithm

    % Number of points
    [~,n1] = size(pt1);
    [~,n2] = size(pt2);
    n = min(n1,n2);

    % Attempt to normalise each set of points so that the origin 
    % is at centroid and mean distance from origin is sqrt(2).
    [pt1, T1] = normalise2dpts(pt1);
    [pt2, T2] = normalise2dpts(pt2);

    % Construct A
    A = zeros(2*n,9);

    for i=1:n

        xip = pt2(1,i);
        yip = pt2(2,i);
        wip = pt2(3,i);    
        point1 = pt1(:,i)';

        A(2*i-1,1:3) = [0 0 0];    %a1
        A(2*i-1,4:6) = -wip*point1;%a2
        A(2*i-1,7:9) = yip*point1; %a3
        A(2*i,1:3) = wip*point1;   %a4
        A(2*i,4:6) = [0 0 0];      %a5
        A(2*i,7:9) = -xip*point1;  %a6
    end

    [U,D,V] = svd(A);
    h = V(:,end);

    H = reshape(h,3,3)';
    H = T2\H*T1;
    H = H / H(3,3);

end


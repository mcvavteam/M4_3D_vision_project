
% Compute the angle between two lines in homogeneous coordinates:
function angle = compute_angle(l1, l2)
    % Cartesian coordinates of the lines:
    l1c = [l1(1) / l1(3), l1(2) / l1(3)]';
    l2c = [l2(1) / l2(3), l2(2) / l2(3)]';
    angle = acos((l1c'*l2c)/sqrt((l1c'*l1c)*(l2c'*l2c)));
    % Pass to degrees:
    angle = angle / (2 * pi) * 360;
    return

end
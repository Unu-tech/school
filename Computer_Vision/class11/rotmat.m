function R = rotmat(ang)
%% 3D rotational matrix
% ang is a 3-component vector contains angles (in radian) about X, Y, Z axes.
Rx = [
    1, 0, 0 ;
    0, cos(ang(1)), -sin(ang(1)) ;
    0, sin(ang(1)),  cos(ang(1)) ] ;
Ry = [
    cos(ang(2)), 0, sin(ang(2)) ;
    0, 1, 0 ;
    -sin(ang(2)), 0, cos(ang(2)) ] ;
Rz = [
    cos(ang(3)), -sin(ang(3)), 0 ;
    sin(ang(3)), cos(ang(3)), 0 ;
    0, 0, 1 ] ;
R = Rz * Ry * Rx ;

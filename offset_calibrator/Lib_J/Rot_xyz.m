function [R] = Rot_xyz(xd,yd,zd)
%ROT_EU 이 함수의 요약 설명 위치
%   자세한 설명 위치
    R = Rotd_axis(1,xd)*Rotd_axis(2,yd)*Rotd_axis(3,zd);
end


function [R] = Rot_xyz(xd,yd,zd)
%ROT_EU �� �Լ��� ��� ���� ��ġ
%   �ڼ��� ���� ��ġ
    R = Rotd_axis(1,xd)*Rotd_axis(2,yd)*Rotd_axis(3,zd);
end


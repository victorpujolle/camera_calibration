function [R] = Rot_zyx(zd,yd,xd)
%ROT_EU �� �Լ��� ��� ���� ��ġ
%   �ڼ��� ���� ��ġ
    R = Rotd_axis(3,zd)*Rotd_axis(2,yd)*Rotd_axis(1,xd);
end


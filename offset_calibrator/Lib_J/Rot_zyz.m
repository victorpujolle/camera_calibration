function [R] = Rot_eu(thetad,phid,psid)
%ROT_EU �� �Լ��� ��� ���� ��ġ
%   �ڼ��� ���� ��ġ
    R = Rotd_axis(3,thetad)*Rotd_axis(2,phid)*Rotd_axis(3,psid);
end


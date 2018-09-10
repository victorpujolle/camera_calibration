function [R] = Rot_eu(thetad,phid,psid)
%ROT_EU 이 함수의 요약 설명 위치
%   자세한 설명 위치
    R = Rotd_axis(3,thetad)*Rotd_axis(2,phid)*Rotd_axis(3,psid);
end


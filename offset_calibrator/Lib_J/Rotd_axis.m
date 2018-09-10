function [ R ] = Rotd_axis( axis, q )
%ROT_AXIS 이 함수의 요약 설명 위치
%   자세한 설명 위치
    switch axis
        case 1
            R = [1,0,0;0,cosd(q),-sind(q);0,sind(q),cosd(q)];
        case 2
            R = [cosd(q),0,sind(q);0,1,0;-sind(q),0,cosd(q)];
        case 3
            R = [cosd(q),-sind(q),0;sind(q),cosd(q),0;0,0,1];
end


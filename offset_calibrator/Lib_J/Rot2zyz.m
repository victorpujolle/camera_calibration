function [z1,y2,z3] = Rot2zyz(R)
%ROT2XYZ 이 함수의 요약 설명 위치
%   자세한 설명 위치sy = sqrt(R(1,1) * R(1,1) +  R(2,1) * R(2,1))
    sy = sqrt(R(1,3)^2 +  R(2,3)^2);
 
    if  sy > 1e-6
        z1 = atan2d(R(2,3) , R(1,3));
        y2 = atan2d(sy,R(3,3));
        z3 = atan2d(R(3,2), -R(3,1));
    else
        z1 = 0;
        y2 = atan2d(sy,R(3,3));
        z3 = atan2d(R(2,1),R(2,2));
    end
end


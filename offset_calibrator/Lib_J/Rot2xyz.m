function [x,y,z] = Rot2xyz(R)
%ROT2XYZ 이 함수의 요약 설명 위치
%   자세한 설명 위치
    sy = sqrt(R(2,3)^2 +  R(3,3)^2);
     
    singular = sy < 1e-6;
 
    if  ~singular
        x = atan2d(-R(2,3) , R(3,3));
        y = atan2d(R(1,3), sy);
        z = atan2d(-R(1,2), R(1,1));
    else
        x = 0;
        y = atan2d(R(1,3), sy);
        z = atan2d(R(3,2), -R(3,1));
    end
end


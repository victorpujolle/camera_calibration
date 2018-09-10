function [z,y,x] = Rot2zyx(R)
%ROT2XYZ �� �Լ��� ��� ���� ��ġ
%   �ڼ��� ���� ��ġsy = sqrt(R(1,1) * R(1,1) +  R(2,1) * R(2,1))
    sy = sqrt(R(1,1) * R(1,1) +  R(2,1) * R(2,1));
 
    if  sy > 1e-6
        x = atan2d(R(3,2) , R(3,3));
        y = atan2d(-R(3,1), sy);
        z = atan2d(R(2,1), R(1,1));
    else
        x = atan2d(-R(2,3), R(2,2));
        y = atan2d(-R(3,1), sy);
        z = 0;   
    end
end


function [points_col, lines] = get_keypoints_object(object)
%GET_KEYPOINTS 이 함수의 요약 설명 위치
%   자세한 설명 위치

[points_row, lines] = get_keypoints_dict(object.dict);

position = object.position;
euler_zyx = object.euler_zyx;
R=Rot_zyx(euler_zyx(1),euler_zyx(2),euler_zyx(3));
N=size(points_row,1);
points_col=(R*points_row');
points_col=points_col+repmat(position',[1,N]);
end
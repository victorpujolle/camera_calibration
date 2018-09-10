function [points_row, lines] = get_keypoints_dict(dict)
%KEYPOINT_LIBRARY 이 함수의 요약 설명 위치
%   자세한 설명 위치
    points_row = dict.mu';
    lines = dict.lines';
end


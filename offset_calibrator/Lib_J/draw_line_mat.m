function [mat] = draw_line_mat(mat,pt1,pt2)
%DRAW_LINE_GAUSSIAN 이 함수의 요약 설명 위치
%   자세한 설명 위치
    x = [pt1(1) pt2(1)];  % x coordinates (running along matrix columns)
    y = [pt1(2) pt2(2)];   % y coordinates (running along matrix rows)
    nPoints = max(abs(diff(x)), abs(diff(y)))+1;  % Number of points in line
    rIndex = round(linspace(y(1), y(2), nPoints));  % Row indices
    cIndex = round(linspace(x(1), x(2), nPoints));  % Column indices
    index = sub2ind(size(mat), rIndex, cIndex);     % Linear indices
    mat(index) = 255;  % Set the line pixels to the max value of 255 for uint8 types

end
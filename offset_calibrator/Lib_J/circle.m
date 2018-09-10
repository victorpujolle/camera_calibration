function [outputArg1,outputArg2] = circle(points_col,r,color,color_txt)
%CIRCLE �� �Լ��� ��� ���� ��ġ
%   �ڼ��� ���� ��ġ
    N=size(points_col,2);
    for i=1:N
        pos=[points_col(:,i)-[r;r];[r*2;r*2]];
        rectangle('Position',pos,'Curvature',[1 1],'EdgeColor',color,'linewidth',2);
        text(pos(1),pos(2),num2str(i),'color',color_txt)
    end
end


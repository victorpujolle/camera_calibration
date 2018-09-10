function [points_px, projection] = project_points_cam(points_col, cam)
%PROJECT_POINTS 이 함수의 요약 설명 위치
%   자세한 설명 위치
    cam_focus=cam.focus;
    cam_dist=cam.dist;
    cam_euler=cam.euler;
    lens_f=cam.lens_f/cam.sensor_dim(1)*cam.resolution(1);
    offset_x=cam.offset_px(1);
    offset_y=cam.offset_px(2);
    
    N=size(points_col,2);
    Teye=eye(4,4);
    Tfocus=Teye;
    Tfocus(1:3,4)=cam_focus;
    Teuler=Teye;
    Teuler(1:3,1:3)=Rot_zxz(cam_euler(1),cam_euler(2),cam_euler(3));
    Tdist=Teye;
    Tdist(3,4)=cam_dist;
    points_se3=[points_col;ones(1,N)];
    Tbc=Tfocus*Teuler*Tdist;
    projection=inv(Tbc)*points_se3;
    if(projection(3,1)<0)
        projection(2,:)=-projection(2,:);
        projection(3,:)=-projection(3,:);
    end
    points_px=[projection(1,:);projection(2,:)]./repmat(projection(3,:),[2,1])*lens_f+repmat([320+offset_x;240+offset_y],[1,N]);
end
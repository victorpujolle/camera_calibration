function [nbytes] = write_pose(cam, ptarfile)
%WRITE_POSE 이 함수의 요약 설명 위치
%   자세한 설명 위치
        M_pose(1:3) = cam.focus;
        M_pose(4) = cam.dist;
        M_pose(5:7) = cam.euler/180*pi;
        M_pose(8) = cam.lens_f;
        M_pose(9:10) = cam.sensor_dim;
        M_pose(11:12) = cam.offset_px;
        pfileID = fopen(ptarfile,'w');
        for i_str=1:size(M_pose,2)
            if i_str==1
                str_p = num2str(M_pose(i_str));
            else
                str_p = [str_p ', ' num2str(M_pose(i_str))];
            end
        end
        nbytes = fprintf(pfileID,[str_p '\r\n']);
        fclose(pfileID);
end


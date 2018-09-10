function [nbytes] = write_config(objects,ctarfile)
%WRITE_CONFIG 이 함수의 요약 설명 위치
%   자세한 설명 위치
    cfileID = fopen(ctarfile,'w');
    for i_obj=1:length(objects)
        str_p = [num2str(objects(i_obj).cls)];
        str_p = [str_p ', ' num2str(objects(i_obj).scale(1)) ', ' num2str(objects(i_obj).scale(2)) ', ' num2str(objects(i_obj).scale(3))];
        str_p = [str_p ', ' num2str(objects(i_obj).position(1)) ', ' num2str(objects(i_obj).position(2)) ', ' num2str(objects(i_obj).position(3))];
        str_p = [str_p ', ' num2str(objects(i_obj).euler_xyz(1)) ', ' num2str(objects(i_obj).euler_xyz(2)) ', ' num2str(objects(i_obj).euler_xyz(3))];
        nbytes = fprintf(cfileID,[str_p '\r\n']);
    end
    fclose(cfileID);
end


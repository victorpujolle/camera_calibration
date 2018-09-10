clear all;

source_path='train_bak/blender3';
config_path = [source_path '/config'];
pose_path = [source_path '/poses'];

cfiles=dir(config_path);
len_ciles=length(cfiles);
for i_c = 1:len_ciles
    filename=cfiles(i_c).name;
    if (contains(filename,'.csv'))
        filename=filename(1:end-4);

        cfilepath=[config_path '/' filename '.csv'];
        M_conf=csvread(cfilepath);
        num_obj=size(M_conf,1);
        cls=zeros(1,num_obj);
        scale=zeros(3,num_obj);
        position=zeros(3,num_obj);
        euler_xyz=zeros(3,num_obj);
        lens_f=zeros(1,num_obj);
        sensor_dim=zeros(2,num_obj);
        offset_px=zeros(2,num_obj);
        for i_obj=1:num_obj
            obj_config = M_conf(i_obj,:);
            cls(i_obj)=obj_config(1);
            scale(:,i_obj)=obj_config(2:4);
            position(:,i_obj)=obj_config(5:7);
            euler_xyz(:,i_obj)=obj_config(8:10);
            lens_f_ref=obj_config(11);
            sensor_dim_ref=obj_config(12:13);
            offset_px_ref=obj_config(14:15);
        end


        pfile=[pose_path '/' filename '.csv'];
        M_pose=csvread(pfile);
        M_pose(8)=lens_f_ref;
        M_pose(9:10)=sensor_dim_ref;
        M_pose(11:12)=offset_px_ref;

        ptarfile=[pose_path '/' filename '.csv'];
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

        ctarfile=[config_path '/' filename '.csv'];
        cfileID = fopen(ctarfile,'w');
        for i_obj=1:num_obj
            str_p = [num2str(cls(i_obj))];
            str_p = [str_p ', ' num2str(scale(1,i_obj)) ', ' num2str(scale(2,i_obj)) ', ' num2str(scale(3,i_obj))];
            str_p = [str_p ', ' num2str(position(1,i_obj)) ', ' num2str(position(2,i_obj)) ', ' num2str(position(3,i_obj))];
            str_p = [str_p ', ' num2str(euler_xyz(1,i_obj)) ', ' num2str(euler_xyz(2,i_obj)) ', ' num2str(euler_xyz(3,i_obj))];
            nbytes = fprintf(cfileID,[str_p '\r\n']);
        end
        fclose(cfileID);
    end
end
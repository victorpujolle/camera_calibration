clear all;

%% First move images, masks and posese to _bak folders.
%% ADD focal	sensorW	sensorH	dpx	dpy
add_cam_setting = true;
lens_f_ref=319.4593/640*3.6;
sensor_dim_ref=[3.6; 2.4];
offset_px_ref=[-15.013; 64.8108];
%% Convert eu_cam from deg to rad
eu_cam_deg2rad = true;
%% Convert bmp to jpg
convert_extension = true;
imfrom='.bmp';
imto='.jpg';
imto_typeonly='jpg';
%% Index mask files
convert_mask_idx = true;
%% Put on floor (z-offset)
set_z_offset = true;
z_off = 0.0055;


source_path='F:/hourglasstensorlfow-master-aug/datasets/config/';
source_files = ["val_01_00", "val_01_01"];
dest_path='F:/hourglasstensorlfow-master-aug/datasets/val_bak/real_';

for i_sf=1:length(source_files)
    source_file=char(source_files(i_sf));
    im_folder = [dest_path source_file];
    ma_folder = [im_folder '/' 'masks'];
    im_bak_folder = [im_folder '/' 'img_bak'];
    ma_bak_folder = [im_folder '/' 'masks_bak'];
    conf_dest_folder = [im_folder '/' 'config'];
    pose_dest_folder = [im_folder '/' 'poses'];
    pose_path=[im_folder '/poses_bak'];
    if(exist(im_bak_folder,'dir')~=7)
        mkdir(im_bak_folder);
    end
    if(exist(ma_bak_folder,'dir')~=7)
        movefile(ma_folder,ma_bak_folder);
        mkdir(ma_folder);
    end
    if(exist(conf_dest_folder,'dir')~=7)
        mkdir(conf_dest_folder);
    end
    if(exist(pose_dest_folder,'dir')~=7)
        mkdir(pose_dest_folder);
    else
        if(exist(pose_path,'dir')~=7)
            movefile(pose_dest_folder,pose_path);
            mkdir(pose_dest_folder);
        end
    end
    cfilepath=[source_path source_file '.csv'];
    objects=read_config(cfilepath);
    for i_obj=1:length(objects)
        if (set_z_offset)
            objects(i_obj).position=put_on_floor(objects(i_obj).scale,objects(i_obj).position,objects(i_obj).euler_xyz/pi*180,z_off);
        end
    end
    
    
%     M_conf=csvread(cfilepath);
% 	num_obj=size(M_conf,1);
%     cls=zeros(1,num_obj);
%     scale=zeros(3,num_obj);
%     position=zeros(3,num_obj);
%     euler_xyz=zeros(3,num_obj);
%     lens_f=zeros(1,num_obj);
%     sensor_dim=zeros(2,num_obj);
%     offset_px=zeros(2,num_obj);
%     for i_obj=1:num_obj
%         obj_config = M_conf(i_obj,:);
%         cls(i_obj)=obj_config(1);
%         scale(:,i_obj)=obj_config(2:4);
%         position(:,i_obj)=obj_config(5:7);
%         euler_xyz(:,i_obj)=obj_config(8:10);
%         if (set_z_offset)
%             position(:,i_obj)=put_on_floor(scale(:,i_obj),position(:,i_obj),euler_xyz(:,i_obj)/pi*180,z_off);
%         end
%     end
    
    im_files=dir(im_folder);
    len_im_files=length(im_files);
    for i_im = 1:len_im_files
        filename=im_files(i_im).name;
        if(convert_extension)
            if (contains(filename,imfrom))
                BMP = imread([im_folder '/' filename]);
                movefile([im_folder '/' filename],[im_bak_folder '/' filename]);
                filename(end-3:end)=imto;
                imwrite(BMP,[im_folder '/' filename],imto_typeonly);   
            end
        end
        if (contains(filename,imto))
            filename=filename(1:end-4);
            pfilepath=[pose_path '/' filename '.csv'];
            M_pose=csvread(pfilepath);
            if(eu_cam_deg2rad)
                M_pose(5:7) = M_pose(5:7)/180*pi;
            end
            cam.focus=M_pose(1:3);
            cam.dist=M_pose(4);
            cam.euler=M_pose(5:7)/pi*180;
            if(add_cam_setting)
                M_pose(8)=lens_f_ref;
                M_pose(9:10)=sensor_dim_ref;
                M_pose(11:12)=offset_px_ref;
                cam.lens_f=M_pose(8);
                cam.sensor_dim=M_pose(9:10);
                cam.offset_px=M_pose(11:12);
                cam.resolution=[640;480];
            end
            ptarfile=[pose_dest_folder '/' filename '.csv'];
            write_pose(cam, ptarfile);
            
            ctarfile=[conf_dest_folder '/' filename '.csv'];
            write_config(objects,ctarfile);
            
            for i_obj=1:length(objects)
                if(convert_mask_idx)
                    from_ma_name=[filename '_' num2str(objects(i_obj).cls+1,'%02i') imfrom];
                    to_ma_name=[filename '_' num2str(i_obj-1,'%02i') '_' num2str(objects(i_obj).cls,'%02i') imto];
                    BMP = imread([ma_bak_folder '/' from_ma_name ]);
                    imwrite(BMP,[ma_folder '/' to_ma_name],imto_typeonly); 
                end
            end
        end
    end
end
clear all;
import_snakeyaml;
imfrom='.png';
imto='.jpg';
imto_extension_only='jpg';

resolution = [640; 480];
sensor_dim= [3.6; 2.4];
Tunit=1/1000;

idx='13';
cadpath = 'F:\dataset\hinterstoisser\models\cad13.mat';
data_path=['F:/dataset/hinterstoisser/test/' idx '/'];
imfolder='rgb/';
dpfolder='depth/';
target_path=['F:/dataset/hinterstoisser/test/' idx '/hinter' idx '/'];
mafolder='masks/';
posefolder='poses/';
configfolder='config/';
cadfolder='cad/';
mkdir(target_path);
mkdir([target_path mafolder]);
mkdir([target_path posefolder]);
mkdir([target_path configfolder]);
mkdir([target_path cadfolder]);

%%% load cad
%%%%% gen_cad first%%%%%
dat=load(cadpath);
cad=dat.cad;
save([target_path cadfolder 'cad.mat'],'cad');


%%% read yml
info = Read_yml_J([data_path 'info.yml']);
gt = Read_yml_J([data_path 'gt.yml']);
dat_num = length(gt.keyset);
for i_dat=1:dat_num
    key=info.keyset(i_dat);
    cam_K=cell2mat(reshape(info.vals{i_dat}.vals{1},[3,3])');
    depth_scale = info.vals{i_dat}.vals{2};
%     view_level = info.vals{i_dat}.vals{3};
    cam_R_m2c = cell2mat(reshape(gt.vals{i_dat}{1}.vals{1},[3,3])');
    cam_t_m2c = cell2mat(gt.vals{i_dat}{1}.vals{2})*Tunit;
    obj_bb = gt.vals{i_dat}{1}.vals{3};
    obj_id = gt.vals{i_dat}{1}.vals{4};
    
    %%%copy_img
    filename = num2str(double(key),'%04i');
    filename6 = num2str(double(key),'%06i');
    impath=[data_path imfolder filename imfrom];
    dppath=[data_path dpfolder filename imfrom];
    im=imread(impath);
    dp=imread(dppath);
    imwrite(im,[target_path filename6 imto],imto_extension_only);
    
    %%%draw_mask
    obj_count=0;
    ma=uint8(dp>0)*255;
    imwrite(ma,[target_path mafolder filename6 '_' num2str(obj_count,'%02i') '_' num2str(obj_id,'%02i') imto],imto_extension_only);
    Tco=[ cam_R_m2c, cam_t_m2c; 0,0,0,1];
    
    
    %%%write_pose
    cam.focus = [0;0;0];
    cam.dist= 0;
    cam.euler= [0;0;0];
    cam.lens_f = mean([cam_K(1,1),cam_K(2,2)])/resolution(1)*sensor_dim(1);
    cam.sensor_dim = sensor_dim;
    cam.offset_px = cam_K(1:2,3)-resolution/2;
    cam.resolution = resolution;
    ptarfile=[target_path posefolder filename6 '.csv'];
    write_pose(cam, ptarfile);
    
    %%%write_config
    [x_,y_,z_]=Rot2xyz(cam_R_m2c);
    object.cls=obj_id;
    object.scale=[0.1,0.1,0.1];
    object.position=cam_t_m2c;
    object.euler_xyz=[x_,y_,z_]/180*pi;
    objects = [object];
    ctarfile=[target_path configfolder filename6 '.csv'];
    write_config(objects,ctarfile);
end
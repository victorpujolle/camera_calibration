clear all;
addpath('./Lib_J')
global calib_data; %% calibration data
global range_dm; %% deviation scale of mask position
global range_dc; %% deviation scale of camera mount offset
global range_df; %% deviation scale of focal length value
global range_dpx; %% deviation scale of pixel offset
global range_angle; %% deviation scale of pixel angular offsets
global marker_rotation; %% marker rotation - marker may positioned with 0, 90, 180, 270 degrees rotation
global T_cb px M0 f;    %% values derived from transform matrix 
                     %% base-to-cam transfrom matrix, detected marker corners
                     %% default marker position, default focal length
                     
savefile='../Offset_list/offset.csv';
calib_data=load('../Calib_list/1.csv');
range_dm = 0.10;
range_dc = 0.05;
range_df = 30;
range_dpx = 30;
range_angle = 10;

minval=100; %initial minimum value
repeat_no=50;
options = optimset('MaxFunEvals',1000000); %minimizer option, iteration number

i_del=boolean(ones(size(calib_data,1),1));
for i=1:size(calib_data,1)
    if((calib_data(i,1)==0) && (calib_data(i,2)==0))
        i_del(i)=true; %calibration data rows without marker detections
    else
        i_del(i)=false;
    end
end
calib_data(i_del,:)=[]; %deleting calibration data rows without detections

process_calib_data; %calculate T_cb px M0 f
for i=1:repeat_no %repeat minimization with many random initialization
    X0=rand(1,16)-0.5; % random initialization (normalized to -0.5 to 0.5
    for j=3 %marker rotation inde
        marker_rotation = 90*(j-1); %j=3: assum marker is 180d rotated 
        X=fminsearch(@calc_std_offset,X0,options); %minimization
        val = calc_std_offset(X); % re-caculate minimum value
        if val<minval
            minval=val;
            minX=X;
            minJ=j;
        end
        disp([num2str(i) '/' num2str(repeat_no) '-' num2str(j) ': val / minval = ' num2str(val) '/' num2str(minval)]);
    end
end
X=minX;
[u, v, mean_u, mean_v, std_u, std_v, mean_err] = calc_offset(X,true,true);



%%%%%save calbrated offsets%%%%%%%%%%%
saveX=zeros(1,19);
%% marker default position
saveX(1:3)=M0(1,:);

%% calibrated focal length
saveX(4:5) = squeeze(f(1,:))+minX(1:2)*range_df;
%% calibrated center pixel (pixel offset)
saveX(6:7) = minX(3:4)*range_dpx;

%%camera offsets
dxc = X(5)*range_dc;
dyc = X(6)*range_dc;
dzc = X(7)*range_dc;
%%camera angle offsets
dac = X(8)*range_angle;
dbc = X(9)*range_angle;
dcc = X(10)*range_angle;

%base offset ( marker offset)
dxm = X(11)*range_dm;
dym = X(12)*range_dm;
dzm = X(13)*range_dm;
%base angle offset ( marker angle offset)
dam = X(14)*range_angle;
dbm = X(15)*range_angle;
dcm = X(16)*range_angle;

%%%%%%%% invese T_cc and T_bb, to compensate it in robot motion %%%%%%
T_cc_i=[Rot_zxz(dac,dbc,dcc), [dxc;dyc;dzc];0,0,0,1];
T_cc = inv(T_cc_i);
[z1,x2,z3]=Rot2zxz(T_cc(1:3,1:3));
saveX(8:10) = T_cc(1:3,4);
saveX(11:13) = [z1,x2,z3];

T_bb_i=[Rot_zyx(dam,dbm,dcm), [dxm;dym;dzm];0,0,0,1];
T_bb = inv(T_bb_i);
[z,y,x]=Rot2zyx(T_cc(1:3,1:3));
saveX(14:16) = T_bb(1:3,4);
saveX(17:19) = [z,y,x];
csvwrite(savefile,saveX);

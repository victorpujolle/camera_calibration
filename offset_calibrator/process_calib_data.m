len_dat=length(calib_data)';
px=reshape(calib_data(:,1:8),[len_dat,2,4]);
px=permute(px,[1,3,2]);
M0=calib_data(:,9:11);
f=calib_data(:,12:13);
F0=calib_data(:,14:16);
eu=calib_data(:,17:19);
dist=calib_data(:,20);
for i_d=1:len_dat
    T_bf(i_d,:,:)=eye(4,4);
    T_ff(i_d,:,:)=eye(4,4);
    T_fc(i_d,:,:)=eye(4,4);

    T_bf(i_d,1:3,4)=F0(i_d,:);
    T_ff(i_d,1:3,1:3)=Rot_zxz(eu(i_d,1),eu(i_d,2),eu(i_d,3));
    T_fc(i_d,1:3,4)=[0;0;dist(i_d)];
    T_bc(i_d,:,:)=squeeze(T_bf(i_d,:,:))*squeeze(T_ff(i_d,:,:))*squeeze(T_fc(i_d,:,:));
    Rt=squeeze(T_bc(i_d,1:3,1:3));
    P=squeeze(T_bc(i_d,1:3,4))';
    T_cb(i_d,:,:)=eye(4,4);
    T_cb(i_d,1:3,1:3)=Rt';
    T_cb(i_d,1:3,4)=-Rt'*P;
end
function std_uv = calc_std_offset(X)
%CALC_STD_OFFSET �� �Լ��� ��� ���� ��ġ
%   �ڼ��� ���� ��ġ
% X = dx,dy,dz,df,dd,theta,phi,psi
    global calib_data;
    global range_dm;
    global range_dc;
    global range_df;
    global range_dpx;
    global range_angle;
    global marker_rotation;
    global T_cb px M0 f;

    %%% multiply normalized calibration variables with scales
    dfx = X(1)*range_df; %focal length
    dfy = X(2)*range_df;
    dcx = X(3)*range_dpx; %pixel offset
    dcy = X(4)*range_dpx;
    dxc = X(5)*range_dc; % camera offset
    dyc = X(6)*range_dc;
    dzc = X(7)*range_dc;
    dac = X(8)*range_angle; % camera angle offset
    dbc = X(9)*range_angle;
    dcc = X(10)*range_angle;
    dxm = X(11)*range_dm; % marker offset
    dym = X(12)*range_dm;
    dzm = X(13)*range_dm;
    dam = X(14)*range_angle; % marker angle offset
    dbm = X(15)*range_angle;
    dcm = X(16)*range_angle;

    ploton=false;
    plot_rotcheck_on=false;
    
    len_dat=length(calib_data)';
    
    
    R_m=Rot_zxz(0,0,marker_rotation);% marker rotation
    T_cc=[Rot_zxz(dac,dbc,dcc), [dxc;dyc;dzc];0,0,0,1];% cam offset
    T_bb=[Rot_zyx(dam,dbm,dcm), [dxm;dym;dzm];0,0,0,1];% base offset
    % corners: marker is 10cm square shape
    corners=zeros(4,3);
    corners(1,1:2)=[-1,1]*0.05;
    corners(2,1:2)=[1,1]*0.05;
    corners(3,1:2)=[1,-1]*0.05;
    corners(4,1:2)=[-1,-1]*0.05;
    
    % camera matrix
    K_c=[f(1)+dfx,0,320+dcx,0;...
        0,f(2)+dfy,240+dcy,0;
        0,0,1,0];
    for I=1:len_dat % data index
        for J=1:4 % marker corner index
            % calculate corner position in cam coordinate
            P_cc(I,J,:)=T_cc*squeeze(T_cb(I,:,:))*T_bb*[(M0(I,:)'+R_m*corners(J,:)');1];
            if(P_cc(I,J,3)<0) % flip if coordinate is inverted (z is negative)
                P_cc(I,J,2)=-P_cc(I,J,2);
                P_cc(I,J,3)=-P_cc(I,J,3);
            end
            
            % calculate corner position in pixels
            x_c(I,J,:)=K_c*squeeze(P_cc(I,J,:));
            px_c(I,J,:)=[x_c(I,J,1),x_c(I,J,2)]/x_c(I,J,3);
        end
    end
    
	% deviation between calculated positions and detected positions
    uv=px_c-px;
	% mean deviation: this is camera pixel offset
    u=reshape(uv(:,:,1),[len_dat*4,1]);
    v=reshape(uv(:,:,2),[len_dat*4,1]);
    
	% standard deviation: of pixel offsets <- this is minimization target
    std_uv=mean(sqrt((u-mean(u)).^2+(v-mean(v)).^2));
    
    if plot_rotcheck_on
        figure(10);
        clf;
        hold on;
        plot([x_o,y_o,z_o]);
        plot(squeeze(P_co_rec(:,1:3,1)),'--','linewidth',2);
        drawnow;
    end
    if ploton
        figure(1);
        clf;
        hold on;
        plot(u);
        plot(v);
        grid on;
        drawnow;
    end
end


function [u, v, mean_u, mean_v, std_u, std_v, mean_err] = calc_offset(X, ploton, plot_rotcheck_on)
%CALC_OFFSET 이 함수의 요약 설명 위치
%   자세한 설명 위치
% X = dx,dy,dz,df,dd,theta,phi,psi
    global calib_data;
    global range_dm;
    global range_dc;
    global range_df;
    global range_dpx;
    global range_angle;
    global marker_rotation;
    global T_cb px M0 f;

    dfx = X(1)*range_df;
    dfy = X(2)*range_df;
    dcx = X(3)*range_dpx;
    dcy = X(4)*range_dpx;
    dxc = X(5)*range_dc;
    dyc = X(6)*range_dc;
    dzc = X(7)*range_dc;
    dac = X(8)*range_angle;
    dbc = X(9)*range_angle;
    dcc = X(10)*range_angle;
    dxm = X(11)*range_dm;
    dym = X(12)*range_dm;
    dzm = X(13)*range_dm;
    dam = X(14)*range_angle;
    dbm = X(15)*range_angle;
    dcm = X(16)*range_angle;

    ploton=false;
    plot_rotcheck_on=false;
    
    len_dat=length(calib_data)';
    
    R_m=Rot_zxz(0,0,marker_rotation);
    T_cc=[Rot_zxz(dac,dbc,dcc), [dxc;dyc;dzc];0,0,0,1];
    T_bb=[Rot_zyx(dam,dbm,dcm), [dxm;dym;dzm];0,0,0,1];
    corners=zeros(4,3);
    corners(1,1:2)=[-1,1]*0.05;
    corners(2,1:2)=[1,1]*0.05;
    corners(3,1:2)=[1,-1]*0.05;
    corners(4,1:2)=[-1,-1]*0.05;
    K_c=[f(1)+dfx,0,320+dcx,0;...
        0,f(2)+dfy,240+dcy,0;
        0,0,1,0];
    for I=1:len_dat
        for J=1:4
            P_cc_rec(I,J,:)=T_cc*squeeze(T_cb(I,:,:))*T_bb*[(M0(I,:)'+R_m*corners(J,:)');1];
            if(P_cc_rec(I,J,3)<0)
                P_cc_rec(I,J,2)=-P_cc_rec(I,J,2);
                P_cc_rec(I,J,3)=-P_cc_rec(I,J,3);
            end
            x_c(I,J,:)=K_c*squeeze(P_cc_rec(I,J,:));
            px_c(I,J,:)=[x_c(I,J,1),x_c(I,J,2)]/x_c(I,J,3);
        end
    end
    uv=px_c-px;
    u=reshape(uv(:,:,1),[len_dat*4,1]);
    v=reshape(uv(:,:,2),[len_dat*4,1]);

    mean_u=mean(u);
    std_u=std(u);
    mean_v=mean(v);
    std_v=std(v);
    mean_err=mean(sqrt((u-mean(u)).^2+(v-mean(v)).^2))

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

        figure(100);
        clf;
        hold on;
%         x_o=calib_data(:,3);
%         y_o=calib_data(:,4);
        plot3([px;px2],[py;py2],u,'o');
        plot3([px;px2],[py;py2],v,'o');
        xlabel('x');
        ylabel('y');
        zlabel('u,v');
    end

    for iiii=1:len_dat
        figure(33); clf; axis([0,640,0,480]); axis equal;
        circle(squeeze(px(iiii,:,:))',3,'r','k');
        hold on; circle((squeeze(px_c(iiii,:,:))-[mean_u,mean_v])',3,'g','k');
        drawnow;
        mkdir('ims');
        saveas(gcf,['ims/' num2str(iiii,'%04i') '.jpg']);
    end
    disp(['df= ( ' num2str(dfx) ', ' num2str(dfy) ' )']);
    disp(['dcxy= ( ' num2str(dcx) ', ' num2str(dcy) ' )']);
    disp(['dXmask= ( ' num2str(dxm) ', ' num2str(dym) ', ' num2str(dzm) ' )']);
    disp(['dRmask= ( ' num2str(dam) ', ' num2str(dbm) ', ' num2str(dcm) ' )']);
    disp(['dXcam= ( ' num2str(dxc) ', ' num2str(dyc) ', ' num2str(dzc) ' )']);
    disp(['dRcam= ( ' num2str(dac) ', ' num2str(dbc) ', ' num2str(dcc) ' )']);
    disp(['mean u: ' num2str(mean_u) ', std u: ' num2str(std_u)]);
    disp(['mean v: ' num2str(mean_v) ', std v: ' num2str(std_v)]);
    disp(['mean pixel error: ' num2str(mean_err)]);
end
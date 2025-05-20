% plot_tilts.m
%
% generate and save plots of the 3 accelerometer channels for each
% instrument
%

clear; close all

load('../pressure_data/geometry');

POBS_dir=dir('../stitched_data/drift_corrected');
f_list={POBS_dir.name}';
file_check=cellfun(@(v)v(1),f_list);
i_list=find(eq(file_check,'P')); % subselects only files begining with 'P'
POBS_list=f_list(i_list);

% first, sort geometry info by depth
[sd,id]=sort(stadepth,'descend');
sla=stalat(id);
slo=stalon(id);
sname=staname(id);
fname=filename(id);

% map from one sorting to the other
for i=1:length(fname)
    if length(fname{i})>6
        fname{i}='POBS7';
    end
    ii(i)=find(strcmp([fname{i} '.mat'],POBS_list));
end
P_list=POBS_list(ii); % should be depth-sorted

for i=1:length(P_list)
    load(['../stitched_data/drift_corrected/' P_list{i}],'dataf');

    t=dataf.tf;
    x=dataf.xtf; x0=x(1); x=x-x0;
    y=dataf.ytf; y0=y(1); y=y-y0;
    z=dataf.ztf; z0=z(1); z=z-z0;
    T=dataf.Ttf;

    a=sqrt(dataf.xtf.^2+dataf.ytf.^2+dataf.ztf.^2); a0=a(1); a=a-a0;

    figure(6); clf
    subplot(311); hold on
    hx=plot(t,x,'linewidth',1);
    ylabel(['x - ' num2str(x0) ' (m/s^2)'])
    datetick('x',3)
    title(sname{i})
    set(gca,'fontsize',12)
    box on; grid on
    yyaxis right
    plot(t,T,'linewidth',0.5)
    ylabel('T (C)')
    yyaxis left
    ha=plot(t,a,'k','linewidth',1);
    legend([hx ha],'a_x','a_t_o_t_a_l')
    subplot(312); hold on
    plot(t,y,'linewidth',1)
    ylabel(['y - ' num2str(y0) ' (m/s^2)'])
    datetick('x',3)
    set(gca,'fontsize',12)
    box on; grid on
    yyaxis right
    plot(t,T,'linewidth',0.5)
    ylabel('T (C)')
    subplot(313); hold on
    plot(t,z,'linewidth',1)
    ylabel(['z - ' num2str(z0) ' (m/s^2)'])
    datetick('x',3)
    set(gca,'fontsize',12)
    box on; grid on
    yyaxis right
    plot(t,T,'linewidth',0.5)
    ylabel('T (C)')

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../figures/stitching/tilt/' sname{i} '_acceleration'],'-dpng','-r300')

    % WW assumes change in a is drift attributable to vertical (x) channel and
    % subtracts it out. I think this is sound, but I do not understand why
    % removing that drift should make the change in x angle equal to the
    % change in horizontal angle, which is also subject to (unknown) drift?

    theta_x=acos((x+x0)/9.81);
    theta_xc=acos((x+x0-a)/9.81);
    theta_y=asin((y+y0)/9.81);
    theta_z=asin((z+z0)/9.81);
    theta_h=sqrt(theta_y.^2+theta_z.^2);

    figure(7); clf
    subplot(311); hold on
    plot(t,theta_x,'linewidth',1)
    plot(t,theta_xc,'linewidth',1)
    ylabel('\Theta_x (rad)')
    datetick('x',3)
    title(sname{i})
    legend('x','x_c_o_r')
    set(gca,'fontsize',12)
    box on; grid on
    subplot(312); hold on
    plot(t,theta_y,'linewidth',1)
    ylabel('\Theta_y (rad)')
    datetick('x',3)
    set(gca,'fontsize',12)
    box on; grid on
    subplot(313); hold on
    plot(t,theta_z,'linewidth',1)
    ylabel('\Theta_z (m/s^2)')
    datetick('x',3)
    set(gca,'fontsize',12)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../figures/stitching/tilt/' sname{i} '_tiltcomponents'],'-dpng','-r300')

    figure(8); clf; hold on
    plot(t,theta_x-theta_x(1),'k','linewidth',2);
    plot(t,theta_xc-theta_xc(1),'k:','linewidth',2);
    plot(t,theta_h-theta_h(1),'c--','linewidth',2)
    ylabel('\Deltatilt (rad)')
    datetick('x',3)
    title(sname{i})
    legend('Vertical channel','Corrected vertical','Horizontal channels')
    set(gca,'fontsize',12)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../figures/stitching/tilt/' sname{i} '_tilt'],'-dpng','-r300')

    if strcmp(sname{i},'POBS-04') || strcmp(sname{i},'POBS-09') || strcmp(sname{i},'POBS-15')
        disp('Cannot determine platform tilt for this instrument')
    else
        %--- calculate vertical displacement of pressure sensor, assuming WW's geometry
        % tiltmeter is oriented with x up, z pointing away from the upper
        % foot ('A'), and y pointing parallel to the edge between the left ('C')
        % and right ('B') feet, in the direction of C (right handed system)

        %--- note that Spahr's configuration is different!
        % the tiltmeter z axis points towards the 'A' foot and y points
        % parallel to the bottom edge of the platform in the direction of
        % B (still right handed). additionally, his has an external
        % pressure reservoir and outlet that is positioned almost directly
        % overtop foot A

        %% I think this means I can basically neglect the effect of B and C sinking
        %% and only have to worry about them if z tilt suggests A lifting up
        
        % geometry
        lt=0.59; % side length of triangular platform [m]
        ht=lt*sqrt(3)/2; % triangular height of platform [m]
        lp1=0.18; % distance from APG to nearest leg [m]
        lp2=lp1*sin(pi/6); % perpendicular distance from APG to adjacent sides [m]

        % net change in tilt from horizontals
        Y=theta_y(end)-theta_y(1);
        Z=theta_z(end)-theta_z(1);

        % map into tilting on platform's feet
        if Y>0 % implies B subsiding
            % B subsiding --> y tilting upwards
            % tilt in 'B' direction that yields observed tilt in y-direction
            By=Y/cos(pi/6);

            % B subsiding --> z tilting downwards
            % tilt induced in z-direction from above B-direction tilt
            Zb=-By*cos(pi/3); % maintains sign relationship

            % how much tilt required from A to total observed tilt in z-direction?
            Za=Z-Zb; % (Z=Za+Zb)

            % z- and A-directions are antiparallel
            Az=Za; % CHECK SIGN AT THIS POINT -- WHY NO NEGATIVE?

            % convert tilting platform into vertical motion of feet
            Bvert=ht*sin(By);
            Avert=ht*sin(Az);

            % scale motion of feet to motion at APG (similar triangles)
            apg_vertB=lp2/ht*Bvert;
            apg_vertA=lp1/ht*Avert;

            apg_vert=apg_vertA+apg_vertB; % total subsidence of APG [m]

            disp([sname{i} ' subsided at least ' num2str(apg_vert) ' m'])
        else % implies C subsiding
            % C subsiding --> y tilting downwards
            % tilt in 'C' direction that yields observed tilt in y-direction
            Cy=-Y/cos(pi/6); % maintains sign relationship

            % C subsiding --> z tilting downwards
            % tilt induced in z-direction from above C-direction tilt
            Zc=-Cy*cos(pi/3); % maintains sign relationship

            % how much tilt required from A to total observed tilt in z-direction?
            Za=Z-Zc; % (Z=Za+Zc)

            % z- and A-directions are antiparallel
            Az=Za; % CHECK SIGN AT THIS POINT -- WHY NO NEGATIVE?

            if Az<0 % indicates foot A moving up, which we are not allowing
                % so, how much must B tilt instead?
                Zb=Za; % get same amount of tilt on z from B instead of A
                Bz=Zb/cos(pi/3);

                % but this also induces upwards tilt on y
            end

            % convert tilting platform into vertical motion of feet
            Bvert=ht*sin(Cy);
            Avert=ht*sin(Az);

            % scale motion of feet to motion at APG (similar triangles)
            apg_vertB=lp2/ht*Bvert;
            apg_vertA=lp1/ht*Avert;

            apg_vert=apg_vertA+apg_vertB; % total subsidence of APG [m]

            disp([sname{i} ' subsided at least ' num2str(apg_vert) ' m'])
        end
    end
end
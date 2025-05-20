% pressure_stitch.m
%
% Load A-0-A pressure data, remove calibration segments, and see how things
% look before making further adjustments.
%

clear; close all

POBS_dir=dir('../pressure_data_Y2/');
POBS_list={POBS_dir.name}';
file_check=cellfun(@(v)v(1),POBS_list);
i_list=find(eq(file_check,'P')); % subselects only files begining with 'P'
for i=i_list'
    % if ~strcmp(POBS_dir(i).name,'POBS09.mat') && ~strcmp(POBS_dir(i).name,'POBS7_co-located-w-QA15.mat')
    %     continue
    % end
    
    load([POBS_dir(i).folder '/' POBS_dir(i).name],'calInfoAll1')
    if isempty(calInfoAll1)
        continue
    elseif length(calInfoAll1.i0)<=4
        continue
    end
    load([POBS_dir(i).folder '/' POBS_dir(i).name])

    if any(data.lNormState1~=data.lNormState2)
        % remove calibration intervals
        t_part=data.ta(data.lNormState1);
        p1_part=data.a1(data.lNormState1);
        T1_part=data.Ta1(data.lNormState1);

        % start with linear interpolation
        t_stitch=data.ta;
        p1_stitch=interp1(t_part,p1_part,t_stitch);
        T1_stitch=interp1(t_part,T1_part,t_stitch);

        % again, for gauge 2
        t_part=data.ta(data.lNormState2);
        p2_part=data.a2(data.lNormState2);
        T2_part=data.Ta2(data.lNormState2);

        p2_stitch=interp1(t_part,p2_part,t_stitch);
        T2_stitch=interp1(t_part,T2_part,t_stitch);
    else
        % remove calibration intervals
        t_part=data.ta(data.lNormState1);
        p1_part=data.a1(data.lNormState1);
        p2_part=data.a2(data.lNormState1);
        T1_part=data.Ta1(data.lNormState1);
        T2_part=data.Ta2(data.lNormState1);

        % start with linear interpolation
        t_stitch=data.ta;
        p1_stitch=interp1(t_part,p1_part,t_stitch);
        p2_stitch=interp1(t_part,p2_part,t_stitch);
        T1_stitch=interp1(t_part,T1_part,t_stitch);
        T2_stitch=interp1(t_part,T2_part,t_stitch);
    end

    % decimation loop
    td=[];
    p1d=[]; p2d=[]; T1d=[]; T2d=[];
    xtd=[]; ytd=[]; ztd=[]; Ttd=[];
    i1 = 1;
    d2 = floor(t_stitch(1))+1;
    while i1<length(t_stitch)
        i2 = find(t_stitch>=d2,1);
        if isempty(i2)
            break
        end
        [segt,segf,~,~] = downsample_uneven(t_stitch(i1:i2-1),[p1_stitch(i1:i2-1),...
            p2_stitch(i1:i2-1),T1_stitch(i1:i2-1),T2_stitch(i1:i2-1),data.xt(i1:i2-1),...
            data.yt(i1:i2-1),data.zt(i1:i2-1),data.Tt(i1:i2-1)],1/24);
        if length(segt)>24
            keyboard
        end
        td=cat(1,td,segt);
        p1d=cat(1,p1d,segf(:,1));
        p2d=cat(1,p2d,segf(:,2));
        T1d=cat(1,T1d,segf(:,3));
        T2d=cat(1,T2d,segf(:,4));
        xtd=cat(1,xtd,segf(:,5));
        ytd=cat(1,ytd,segf(:,6));
        ztd=cat(1,ztd,segf(:,7));
        Ttd=cat(1,Ttd,segf(:,8));
        i1=i2;
        d2=floor(t_stitch(i2))+1;
    end

    % tidal filter
    tf=td;
    p1f=Z_godin(p1d); p2f=Z_godin(p2d); T1f=Z_godin(T1d); T2f=Z_godin(T2d);
    xtf=Z_godin(xtd); ytf=Z_godin(ytd); ztf=Z_godin(ztd); Ttf=Z_godin(Ttd);

    % remove NaNs from tidal filter
    inan=isnan(p1f);
    tf(inan)=[]; T1f(inan)=[]; T2f(inan)=[]; p2f(inan)=[]; p1f(inan)=[];
    xtf(inan)=[]; ytf(inan)=[]; ztf(inan)=[]; Ttf(inan)=[];

    % interpolate onto monotonic time basis
    % merge into single structure
    dataf=[];
    tff=(tf(1)+datenum(0,0,0,0,30,0):1/24:tf(end)-datenum(0,0,0,0,30,0))';
    dataf.tf=tff; dataf.p1f=interp1(tf,p1f,tff); dataf.p2f=interp1(tf,p2f,tff);
    dataf.T1f=interp1(tf,T1f,tff); dataf.T2f=interp1(tf,T2f,tff);
    dataf.xtf=interp1(tf,xtf,tff); dataf.ytf=interp1(tf,ytf,tff);
    dataf.ztf=interp1(tf,ztf,tff); dataf.Ttf=interp1(tf,Ttf,tff);

    save(['../stitched_data_Y2/' POBS_dir(i).name],'dataf')

    % also save 1 Hz data
    data=[];
    data.t=t_stitch;
    data.p1=p1_stitch;
    data.T1=T1_stitch;
    data.p2=p2_stitch;
    data.T2=T2_stitch;
    
    save(['../stitched_data_Y2/stitched_1Hz/' POBS_dir(i).name(1:end-4)],'data')

    %% PLOTS

    figure(23); clf;
    subplot(211); hold on
    plot(dataf.tf,dataf.p1f-dataf.p1f(end),'linewidth',1)
    barcor1=calInfoAll1.pCal-barInfoAll.pCal;
    plot(barInfoAll.t0p,barcor1-barcor1(end),'o','markersize',10,'linewidth',2)
    ylabel('P (hPa)')
    legend('P','cals','location','northeast')
    datetick('x',6)
    title('Gauge 1')
    set(gca,'fontsize',14)
    box on; grid on
    lim1=ylim;

    subplot(212); hold on
    plot(dataf.tf,dataf.p2f-dataf.p2f(end),'linewidth',1)
    barcor2=calInfoAll2.pCal-barInfoAll.pCal;
    plot(barInfoAll.t0p,barcor2-barcor2(end),'o','markersize',10,'linewidth',2)
    ylabel('P (hPa)')
    legend('P','cals','location','northeast')
    datetick('x',6)
    title('Gauge 2')
    set(gca,'fontsize',14)
    box on; grid on
    lim2=ylim;

    % equilibrate ylimits
    subplot(211)
    ylim([min(lim1(1),lim2(1)) max(lim1(2),lim2(2))])
    subplot(212)
    ylim([min(lim1(1),lim2(1)) max(lim1(2),lim2(2))])

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../figures_Y2/stitching/' POBS_dir(i).name(1:end-4) '_stitched'],'-dpng','-r300')
end
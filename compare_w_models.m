% compare_w_models.m
%
% Compares A-0-A drift-corrected pressure data against oceanographic
% models, other oceanographic observables
%

clear; close all

% load HYCOM pressures
load('../../hikurangi/ocean_data/HYCOM/hycom_pressure.mat')
% load GLORYS pressures
load('../../hikurangi/ocean_data/GLORYS/glorys_pressure_update.mat')
% load satellite altimetry
load('../../hikurangi/ocean_data/SSH/satellite_ssh_update.mat')
% load ECCO2 pressures
load('../../hikurangi/ocean_data/ECCO2/ecco2_pressure_update.mat')

POBS_dir=dir('../stitched_data/drift_corrected');
POBS_list={POBS_dir.name}';
file_check=cellfun(@(v)v(1),POBS_list);
i_list=find(eq(file_check,'P')); % subselects only files begining with 'P'

% i_list=[13 12 11 5 15 9 7 8 14 10]'; % inelegant fix to get them in depth order
i_list=[13 12 11 5 6 15 9 7 8 14 10 4]'; % includes 09 and 10

% low-pass filter parameters
[b,a]=butter(4,(2/(30*24)),'low'); % assumes hourly sample rate

space=5; % vertical spacing for stacked plots [hPa]

h_rate=nan(12,1); hf_rate=nan(12,1);
s_rate=nan(12,1); sf_rate=nan(12,1);
g_rate=nan(12,1); gf_rate=nan(12,1);
a_rate=nan(12,1); af1_rate=nan(12,1); af2_rate=nan(12,1);
e_rate=nan(12,1); ef_rate=nan(12,1);

n=1;
for i=1:length(i_list)
    % parameter for helping stacked plots display better
    if i<=4
        xspace=0;
    else
        xspace=1;
    end

    k=i_list(i);
    load([POBS_dir(k).folder '/' POBS_dir(k).name])

    %--- optionally cut out first bit of data with poor calibrations
    % 2 months cut (Laura's preference)
    % dataf.tf(1:60*24)=[]; dataf.p1_dcor(1:60*24)=[]; dataf.p2_dcor(1:60*24)=[];
    % 2 weeks cut (my preference)
    dataf.tf(1:14*24)=[]; dataf.p1_dcor(1:14*24)=[]; dataf.p2_dcor(1:14*24)=[];

    if strcmp(POBS_dir(k).name,'POBS12.mat') % gauge 1 data were truncated with NaNs, which will cause problems later
        inan=find(isnan(dataf.p1_dcor));
        dataf.p1_dcor(inan)=dataf.p2_dcor(inan)-(dataf.p2_dcor(inan(1))-dataf.p1_dcor(inan(1)-1));
    end

    % identifier for station
    shortname=POBS_dir(k).name(1:end-4);
    if length(shortname)==5
        shortname=[shortname(1:4) '-0' shortname(5)];
    elseif length(shortname)==6
        shortname=[shortname(1:4) '-' shortname(5:6)];
    else
        shortname=[shortname(1:4) '-0' shortname(5)];
    end
    
    % get linear fits to APGs
    pm=polyfit(dataf.tf-dataf.tf(1),dataf.p1_dcor-mean(dataf.p1_dcor),1);
    p1_l=polyval(pm,dataf.tf-dataf.tf(1));
    pm=polyfit(dataf.tf-dataf.tf(1),dataf.p2_dcor-mean(dataf.p2_dcor),1);
    p2_l=polyval(pm,dataf.tf-dataf.tf(1));
    a_rate(i)=pm(1)*365; % store A-0-A rate
    
    % %% HYCOM
    % % find model index that corresponds to current station
    % i_hy=find(strcmp(shortname,p_hy.nearest_sta));
    % 
    % % trim HYCOM data to match APG
    % th=p_hy.t{i_hy};
    % ph=p_hy.bpa{i_hy};
    % ikeep=find(th>=(dataf.tf(1)) & th<=dataf.tf(end));
    % ikeep=[ikeep(1)-1;ikeep;ikeep(end)+1];
    % ph=ph(ikeep);
    % th=th(ikeep);
    % 
    % % linear fits to the time series
    % pm=polyfit(th-th(1),ph,1);
    % h_rate(i)=pm(1)*365; % store HYCOM rate
    % ph_l=polyval(pm,th-th(1));
    % 
    % figure(45); clf; hold on
    % subplot(211); hold on
    % plot(dataf.tf,dataf.p1_dcor-mean(dataf.p1_dcor),'linewidth',1)
    % plot(th,ph,'linewidth',1)
    % plot(dataf.tf,p1_l,'color',[0 114 189]/255,'linewidth',1)
    % plot(th,ph_l,'r','linewidth',1)
    % ylabel('P (hPa)')
    % legend('APG','HYCOM','location','northwest')
    % datetick('x',6)
    % title({shortname; 'Gauge 1'})
    % set(gca,'fontsize',14)
    % box on; grid on
    % subplot(212); hold on
    % plot(dataf.tf,dataf.p2_dcor-mean(dataf.p2_dcor),'linewidth',1)
    % plot(th,ph,'linewidth',1)
    % plot(dataf.tf,p2_l,'color',[0 114 189]/255,'linewidth',1)
    % plot(th,ph_l,'r','linewidth',1)
    % ylabel('P (hPa)')
    % legend('APG','HYCOM','location','northwest')
    % datetick('x',6)
    % title('Gauge 2')
    % set(gca,'fontsize',14)
    % box on; grid on
    % 
    % fh=gcf;
    % fh.PaperUnits='inches';
    % fh.PaperPosition=[0 0 8.5 11];
    % % print(['../figures/trend_comps/hycom/hycom_' POBS_dir(k).name(1:end-4)],'-dpng','-r300')
    % 
    % %----- stacked timeseries of APG-HYCOM comparison
    % ph=interp1(th,ph,dataf.tf);
    % 
    % % apply filter defined at start of file
    % p1f=filtfilt(b,a,dataf.p1_dcor);
    % p2f=filtfilt(b,a,dataf.p2_dcor);
    % phf=filtfilt(b,a,ph);
    % 
    % % optional scaling
    % phf=phf-mean(phf);
    % p1f=p1f-mean(p1f); p2f=p2f-mean(p2f);
    % mh=inv(phf'*phf)*phf'*p2f;
    % % disp(['HYCOM:' newline shortname newline num2str(mh) newline])
    % disp([shortname newline 'HYCOM: '  num2str(mh)])
    % 
    % mh=1; % comment out to use least-squares scaling
    % 
    % % find new linear fits
    % pm=polyfit(dataf.tf-dataf.tf(1),phf*mh,1);
    % hf_rate(i)=pm(1)*365; % store HYCOM rate
    % pm=polyfit(dataf.tf-dataf.tf(1),p1f,1);
    % af1_rate(i)=pm(1)*365; % store Gauge 1 rate
    % pm=polyfit(dataf.tf-dataf.tf(1),p2f,1);
    % af2_rate(i)=pm(1)*365; % store Gauge 2 rate
    % 
    % pcor=p2f-phf*mh;
    % figure(114); hold on
    % plot(dataf.tf,p2f+space*n+xspace,'color',[0 114 189]/255,'linewidth',1)
    % plot(dataf.tf,phf*mh+space*n+xspace,'r','linewidth',1)
    % plot(dataf.tf,pcor+space*n+xspace,'k','linewidth',2)
    % text(dataf.tf(1)-60,space*n+xspace,{shortname;[num2str(round(af2_rate(i),1)) ' hPa/yr'];...
    %     [num2str(round(af2_rate(i)-hf_rate(i),1)) ' hPa/yr']},'fontsize',12)

    % %% SATELLITE SSH
    % % find grid index that corresponds to current station
    % i_ssh=find(strcmp(shortname,sat_ssh.nearest_sta));
    % 
    % % trim satellite data to match APG
    % ts=sat_ssh.t{i_ssh};
    % ps=sat_ssh.sla_int{i_ssh};
    % ikeep=find(ts>=(dataf.tf(1)) & ts<=dataf.tf(end));
    % ikeep=[ikeep(1)-1;ikeep;ikeep(end)+1];
    % ps=ps(ikeep);
    % ts=ts(ikeep);
    % 
    % % linear fits to the time series
    % pm=polyfit(ts-ts(1),ps,1);
    % s_rate(i)=pm(1)*365; % store satellite rate
    % ps_l=polyval(pm,ts-ts(1));
    % 
    % figure(46); clf; hold on
    % subplot(211); hold on
    % plot(dataf.tf,dataf.p1_dcor-mean(dataf.p1_dcor),'linewidth',1)
    % plot(ts,ps,'linewidth',1)
    % plot(dataf.tf,p1_l,'color',[0 114 189]/255,'linewidth',1)
    % plot(ts,ps_l,'r','linewidth',1)
    % ylabel('P (hPa)')
    % legend('APG','SSH','location','northwest')
    % datetick('x',6)
    % title({shortname; 'Gauge 1'})
    % set(gca,'fontsize',14)
    % box on; grid on
    % subplot(212); hold on
    % plot(dataf.tf,dataf.p2_dcor-mean(dataf.p2_dcor),'linewidth',1)
    % plot(ts,ps,'linewidth',1)
    % plot(dataf.tf,p2_l,'color',[0 114 189]/255,'linewidth',1)
    % plot(ts,ps_l,'r','linewidth',1)
    % ylabel('P (hPa)')
    % legend('APG','SSH','location','northwest')
    % datetick('x',6)
    % title('Gauge 2')
    % set(gca,'fontsize',14)
    % box on; grid on
    % 
    % fh=gcf;
    % fh.PaperUnits='inches';
    % fh.PaperPosition=[0 0 8.5 11];
    % % print(['../figures/trend_comps/satellite/satellite_' POBS_dir(k).name(1:end-4)],'-dpng','-r300')
    % 
    % %----- stacked timeseries of APG-SSH comparison
    % ps=interp1(ts,ps,dataf.tf);
    % 
    % % apply filter defined at start of file
    % p1f=filtfilt(b,a,dataf.p1_dcor);
    % p2f=filtfilt(b,a,dataf.p2_dcor);
    % psf=filtfilt(b,a,ps);
    % 
    % % optional scaling
    % psf=psf-mean(psf);
    % p1f=p1f-mean(p1f); p2f=p2f-mean(p2f);
    % ms=inv(psf'*psf)*psf'*p2f;
    % % disp(['HYCOM:' newline shortname newline num2str(mh) newline])
    % disp([shortname newline 'SSH: '  num2str(ms)])
    % 
    % ms=1; % comment out to use least-squares scaling
    % 
    % % find new linear fits
    % pm=polyfit(dataf.tf-dataf.tf(1),psf*ms,1);
    % sf_rate(i)=pm(1)*365; % store SSH rate
    % pm=polyfit(dataf.tf-dataf.tf(1),p1f,1);
    % af1_rate(i)=pm(1)*365; % store Gauge 1 rate
    % pm=polyfit(dataf.tf-dataf.tf(1),p2f,1);
    % af2_rate(i)=pm(1)*365; % store Gauge 2 rate
    % 
    % pcor=p2f-psf*ms;
    % figure(113); hold on
    % plot(dataf.tf,p2f+3*space*n+xspace,'color',[0 114 189]/255,'linewidth',1)
    % plot(dataf.tf,psf*ms+3*space*n+xspace,'r','linewidth',1)
    % plot(dataf.tf,pcor+3*space*n+xspace,'k','linewidth',2)
    % if strcmp(shortname,'POBS-09')
    %     text(dataf.tf(1)-60,3*space*n+xspace,{shortname;[num2str(round(af2_rate(i),1)) ' hPa/yr'];...
    %         [num2str(round(af2_rate(i)-sf_rate(i),1)) ' hPa/yr']},'fontsize',12)
    % else
    %     text(dataf.tf(end)+10,3*space*n+xspace,{shortname;[num2str(round(af2_rate(i),1)) ' hPa/yr'];...
    %         [num2str(round(af2_rate(i)-sf_rate(i),1)) ' hPa/yr']},'fontsize',12)
    % end

    %% GLORYS
    % find model index that corresponds to current station
    i_glo=find(strcmp(shortname,p_glo.nearest_sta));

    % trim GLORYS data to match APG
    tg=p_glo.t{i_glo};
    pg=p_glo.bpa{i_glo};
    ikeep=find(tg>=(dataf.tf(1)) & tg<=dataf.tf(end));
    if ikeep(end)==length(tg)
        ng=ceil(dataf.tf(end)-tg(end));
        ng=(1:ng)';
        pg=interp1(tg,pg,[tg;tg(end)+ng],'linear','extrap');
        tg=[tg;tg(end)+ng];
        ikeep=find(tg>=(dataf.tf(1)) & tg<=dataf.tf(end));
    end
    ikeep=[ikeep(1)-1;ikeep;ikeep(end)+1];
    pg=pg(ikeep);
    tg=tg(ikeep);

    % linear fits to the time series
    pm=polyfit(tg-tg(1),pg,1);
    g_rate(i)=pm(1)*365; % store GLORYS rate
    pg_l=polyval(pm,tg-tg(1));

    figure(47); clf; hold on
    subplot(211); hold on
    plot(dataf.tf,dataf.p1_dcor-mean(dataf.p1_dcor),'linewidth',1)
    plot(tg,pg,'linewidth',1)
    plot(dataf.tf,p1_l,'color',[0 114 189]/255,'linewidth',1)
    plot(tg,pg_l,'r','linewidth',1)
    ylabel('P (hPa)')
    legend('APG','GLORYS','location','northwest')
    datetick('x',6)
    title({shortname; 'Gauge 1'})
    set(gca,'fontsize',14)
    box on; grid on
    subplot(212); hold on
    plot(dataf.tf,dataf.p2_dcor-mean(dataf.p2_dcor),'linewidth',1)
    plot(tg,pg,'linewidth',1)
    plot(dataf.tf,p2_l,'color',[0 114 189]/255,'linewidth',1)
    plot(tg,pg_l,'r','linewidth',1)
    ylabel('P (hPa)')
    legend('APG','GLORYS','location','northwest')
    datetick('x',6)
    title('Gauge 2')
    set(gca,'fontsize',14)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    % print(['../figures/trend_comps/glorys/glorys_' POBS_dir(k).name(1:end-4)],'-dpng','-r300')

    %----- stacked timeseries of APG-GLORYS comparison
    pg=interp1(tg,pg,dataf.tf);

    % apply filter defined at start of file
    p1f=filtfilt(b,a,dataf.p1_dcor);
    p2f=filtfilt(b,a,dataf.p2_dcor);
    pgf=filtfilt(b,a,pg);

    % optional scaling
    pgf=pgf-mean(pgf);
    p1f=p1f-mean(p1f); p2f=p2f-mean(p2f);
    mg=inv(pgf'*pgf)*pgf'*p2f;
    % disp(['HYCOM:' newline shortname newline num2str(mh) newline])
    disp([shortname newline 'GLORYS: '  num2str(mg)])

    mg=1; % comment out to use least-squares scaling

    % find new linear fits
    pm=polyfit(dataf.tf-dataf.tf(1),pgf*mg,1);
    gf_rate(i)=pm(1)*365; % store GLORYS rate
    pm=polyfit(dataf.tf-dataf.tf(1),p1f,1);
    af1_rate(i)=pm(1)*365; % store Gauge 1 rate
    pm=polyfit(dataf.tf-dataf.tf(1),p2f,1);
    af2_rate(i)=pm(1)*365; % store Gauge 2 rate

    pcor=p2f-pgf*mg;
    figure(112); hold on
    plot(dataf.tf,p2f+space*n+xspace,'color',[0 114 189]/255,'linewidth',1)
    plot(dataf.tf,pgf*mg+space*n+xspace,'r','linewidth',1)
    plot(dataf.tf,pcor+space*n+xspace,'k','linewidth',2)
    if strcmp(shortname,'POBS-09')
        text(dataf.tf(1)-60,space*n+xspace,{shortname;[num2str(round(af2_rate(i),1)) ' hPa/yr'];...
            [num2str(round(af2_rate(i)-gf_rate(i),1)) ' hPa/yr']},'fontsize',12)
    else
        text(dataf.tf(end)+10,space*n+xspace,{shortname;[num2str(round(af2_rate(i),1)) ' hPa/yr'];...
            [num2str(round(af2_rate(i)-gf_rate(i),1)) ' hPa/yr']},'fontsize',12)
    end

    %% ECCO2

    % if strcmp(shortname,'POBS-10') || strcmp(shortname,'POBS-09')
    %     n=n+1;
    %     continue
    % end

    % find model index that corresponds to current station
    if strcmp(shortname,'POBS-10')
        i_eco=find(strcmp('POBS-05',ecco2.nearest_sta));
    elseif strcmp(shortname,'POBS-09')
        i_eco=find(strcmp('POBS-16',ecco2.nearest_sta));
    else
        i_eco=find(strcmp(shortname,ecco2.nearest_sta));
    end

    % trim ECCO data to match APG
    te=ecco2.t{i_eco};
    pe=ecco2.bpa_int{i_eco};
    ikeep=find(te>=(dataf.tf(1)) & te<=dataf.tf(end));
    if strcmp(shortname,'POBS-09')
        ikeep=[ikeep(1)-1;ikeep];
    else
        ikeep=[ikeep(1)-1;ikeep;ikeep(end)+1];
    end
    pe=pe(ikeep);
    te=te(ikeep);

    % linear fits to the time series
    pm=polyfit(te-te(1),pe,1);
    e_rate(i)=pm(1)*365; % store ECCO2 rate
    pe_l=polyval(pm,te-te(1));

    figure(49); clf; hold on
    subplot(211); hold on
    plot(dataf.tf,dataf.p1_dcor-mean(dataf.p1_dcor),'linewidth',1)
    plot(te,pe,'linewidth',1)
    plot(dataf.tf,p1_l,'color',[0 114 189]/255,'linewidth',1)
    plot(te,pe_l,'r','linewidth',1)
    ylabel('P (hPa)')
    legend('APG','ECCO2','location','northwest')
    datetick('x',6)
    title({shortname; 'Gauge 1'})
    set(gca,'fontsize',14)
    box on; grid on
    subplot(212); hold on
    plot(dataf.tf,dataf.p2_dcor-mean(dataf.p2_dcor),'linewidth',1)
    plot(te,pe,'linewidth',1)
    plot(dataf.tf,p2_l,'color',[0 114 189]/255,'linewidth',1)
    plot(te,pe_l,'r','linewidth',1)
    ylabel('P (hPa)')
    legend('APG','ECCO2','location','northwest')
    datetick('x',6)
    title('Gauge 2')
    set(gca,'fontsize',14)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    % print(['../figures/trend_comps/ecco2/ecco2_' POBS_dir(k).name(1:end-4)],'-dpng','-r300')

    %----- stacked timeseries of APG-ECCO comparison
    if strcmp(shortname,'POBS-09')
        t1=dataf.tf;
        p1=dataf.p1_dcor;
        p2=dataf.p2_dcor;
        icut=t1>te(end);
        t1=t1(~icut); p1=p1(~icut); p2=p2(~icut); 

        pe=interp1(te,pe,t1);

        % apply filter defined at start of file
        p1f=filtfilt(b,a,dataf.p1_dcor);
        p1f_cut=filtfilt(b,a,p1);
        p2f=filtfilt(b,a,dataf.p2_dcor);
        p2f_cut=filtfilt(b,a,p2);
        pef=filtfilt(b,a,pe);

        % optional scaling
        pef=pef-mean(pef);
        p1f=p1f-mean(p1f); p2f=p2f-mean(p2f);
        p1f_cut=p1f_cut-mean(p1f_cut); p2f_cut=p2f_cut-mean(p2f_cut);
        me=inv(pef'*pef)*pef'*p2f_cut;
        % disp(['HYCOM:' newline shortname newline num2str(mh) newline])
        disp([shortname newline 'ECCO: '  num2str(me)])

        me=1; % comment out to use least-squares scaling

        % find new linear fits
        pm=polyfit(t1-t1(1),pef*me,1);
        ef_rate(i)=pm(1)*365; % store ECCO rate
        pm=polyfit(t1-t1(1),p1f_cut,1);
        af1_rate(i)=pm(1)*365; % store Gauge 1 rate
        pm=polyfit(t1-t1(1),p2f_cut,1);
        af2_rate(i)=pm(1)*365; % store Gauge 2 rate

        pcor=p2f_cut-pef*me;
        figure(111); hold on
        plot(dataf.tf,p2f+space*n+xspace,'color',[0 114 189]/255,'linewidth',1)
        plot(t1,pef*me+space*n+xspace,'r','linewidth',1)
        plot(t1,pcor+space*n+xspace,'k','linewidth',2)
        text(dataf.tf(end)+10,space*n+xspace,{shortname;[num2str(round(af2_rate(i),1)) ' hPa/yr'];...
            [num2str(round(af2_rate(i)-ef_rate(i),1)) ' hPa/yr']},'fontsize',12)
    else
        pe=interp1(te,pe,dataf.tf);

        % apply filter defined at start of file
        p1f=filtfilt(b,a,dataf.p1_dcor);
        p2f=filtfilt(b,a,dataf.p2_dcor);
        pef=filtfilt(b,a,pe);

        % optional scaling
        pef=pef-mean(pef);
        p1f=p1f-mean(p1f); p2f=p2f-mean(p2f);
        me=inv(pef'*pef)*pef'*p2f;
        % disp(['HYCOM:' newline shortname newline num2str(mh) newline])
        disp([shortname newline 'ECCO: '  num2str(me)])

        me=1; % comment out to use least-squares scaling

        % find new linear fits
        pm=polyfit(dataf.tf-dataf.tf(1),pef*me,1);
        ef_rate(i)=pm(1)*365; % store ECCO rate
        pm=polyfit(dataf.tf-dataf.tf(1),p1f,1);
        af1_rate(i)=pm(1)*365; % store Gauge 1 rate
        pm=polyfit(dataf.tf-dataf.tf(1),p2f,1);
        af2_rate(i)=pm(1)*365; % store Gauge 2 rate

        pcor=p2f-pef*me;
        figure(111); hold on
        plot(dataf.tf,p2f+space*n+xspace,'color',[0 114 189]/255,'linewidth',1)
        plot(dataf.tf,pef*me+space*n+xspace,'r','linewidth',1)
        plot(dataf.tf,pcor+space*n+xspace,'k','linewidth',2)
        text(dataf.tf(end)+10,space*n+xspace,{shortname;[num2str(round(af2_rate(i),1)) ' hPa/yr'];...
            [num2str(round(af2_rate(i)-ef_rate(i),1)) ' hPa/yr']},'fontsize',12)
    end

    n=n+1;
end

dirlist={'ecco2','glorys','satellite','hycom'};
for j=1:4
    figure(110+j)
    ylabel('P (hPa)')
    legend('APG',dirlist{j},'Difference','location','northeast')
    xlim([datenum(2022,10,01) datenum(2024,04,15)])
    datetick('x',6,'keeplimits')
    set(gca,'fontsize',14)
    ylim([-5 75])
    set(gca,'ytick',-5:5:75)
    box on; grid on
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../figures/trend_comps/' dirlist{j} '/' dirlist{j} '_filteredstack'],'-dpng','-r300')
    print(['../figures/trend_comps/' dirlist{j} '/' dirlist{j} '_filteredstack'],'-depsc','-vector')
end

%% MAPS OF OCEAN-INDUCED CHANGE

addpath('../../hikurangi/code/m_map')
load('../pressure_data/geometry');

% reconcile different station order
stalist = {POBS_dir(i_list).name};
for i=1:length(stalist)
    staID = str2double(stalist{(i)}(5:end-4));
    if isnan(staID)
        staID=7;
    end
    staID_padded = sprintf('%02d',staID);
    ii(i) = find(strcmp(['POBS-' staID_padded],staname));
end

% setup base map
figure(54); clf; hold on
lat1=-40; lat2=-38.5;
lon1=177.5; lon2=179.5;
m_proj('mercator','longitudes',[lon1 lon2],'latitudes',[lat1 lat2]);

X=linspace(174.5,180,1320);
Y=linspace(-37.5,-42,1080);
Z=imread('../../hikurangi/gshhg-bin-2.3.6/NZ_bathymetry.tiff');
Z(Z>=0)=NaN;
[~,hc]=m_contourf(X,Y,Z,-5000:100:0,'linecolor','none');

m_gshhs_i('patch',[.7 .7 .7]);

m_contour(X,Y,Z,[-150 -150],'k','linewidth',1) % shelf at 150 m?
m_text(177.75,-39.55,'150 m','fontsize',12,'rotation',35)
m_contour(X,Y,Z,[-2750 -2750],'k','linewidth',1) % 2750 m does a decent job visually
m_text(178.4,-39.9,'2750 m','fontsize',12,'rotation',80)

colormap(gca,cmocean('-deep'))
cb1 = colorbar(gca,'eastoutside');
ylabel(cb1,'Depth (m)')
set(cb1,'fontsize',14)
m_grid('xlabeldir','end','fontsize',14);

% station markers
m_plot(stalon,stalat,'^k','markersize',10)
m_text(stalon-0.21,stalat,staname,'fontsize',12)

% scale arrow
m_quiver(177.6,-39.95,0,2/20,'off','color','k','linewidth',2)
m_text(177.65,-39.8,'2 cm/yr','fontsize',12)

%---BASE
h1=m_quiver(stalon(ii),stalat(ii),zeros(12,1),-(a_rate)/20,'off','color','c','linewidth',2);
title('A0A inferred deformation')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/trend_comps/base_def','-dpng','-r300')
print('../figures/trend_comps/base_def','-depsc','-vector')

%---HYCOM
delete(h1)
h1=m_quiver(stalon(ii),stalat(ii),zeros(10,1),-h_rate/20,'off','color','r','linewidth',2);
title('HYCOM inferred deformation')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/trend_comps/hycom/mapview','-dpng','-r300')

%---SATELLITE
delete(h1)
h1=m_quiver(stalon(ii),stalat(ii),zeros(10,1),-s_rate/20,'off','color','r','linewidth',2);
title('Altimetry inferred deformation')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/trend_comps/satellite/mapview','-dpng','-r300')

%---GLORYS
delete(h1)
h1=m_quiver(stalon(ii),stalat(ii),zeros(12,1),-g_rate/20,'off','color','r','linewidth',2);
title('GLORYS inferred deformation')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/trend_comps/glorys/mapview','-dpng','-r300')

% this one actually compares well, so let's take a difference
delete(h1)
h1=m_quiver(stalon(ii),stalat(ii),zeros(12,1),-(a_rate-g_rate)/20,'off','color','r','linewidth',2);
title('A0A - GLORYS inferred deformation')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/trend_comps/glorys/mapview_dif','-dpng','-r300')
print('../figures/trend_comps/glorys/mapview_dif','-depsc','-vector')

%---ECCO2
delete(h1)
h1=m_quiver(stalon(ii),stalat(ii),zeros(10,1),-e_rate/20,'off','color','r','linewidth',2);
title('ECCO2 inferred deformation')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/trend_comps/ecco2/mapview','-dpng','-r300')

% this one actually compares well, so let's take a difference
delete(h1)
h1=m_quiver(stalon(ii),stalat(ii),zeros(12,1),-(a_rate-e_rate)/20,'off','color','r','linewidth',2);
title('A0A - ECCO2 inferred deformation')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/trend_comps/ecco2/mapview_dif','-dpng','-r300')
print('../figures/trend_comps/ecco2/mapview_dif','-depsc','-vector')
% inspect_postcal_transient.m
%
% Align all data segments relative to calibration times and see if stacking
% reveals a repeatable transient.
%

clear; close all

load('../pressure_data/geometry.mat')
path(path,'../../EKF_spotl/matlab')

POBS_dir=dir('../stitched_data/');
POBS_list={POBS_dir.name}';
file_check=cellfun(@(v)v(1),POBS_list);
i_list=find(~eq(file_check,'.'));
i_list(end)=[];

% reconcile different station order
stalist = {POBS_dir.name};
stalist = stalist(i_list);
for i=1:length(stalist)
    staID = str2double(stalist{(i)}(5:end-4));
    if isnan(staID)
        staID=7;
    end
    staID_padded = sprintf('%02d',staID);
    ii(i) = find(strcmp(['POBS-' staID_padded],staname));
end
staname=staname(ii);
stalat=stalat(ii);
stalon=stalon(ii);
stadepth=stadepth(ii);

for k=1:length(i_list)
    i=i_list(k);
    % % UNCOMMENT TO RUN ONLY SPECIFIED STATIONS
    % if ~strcmp(POBS_dir(i).name,'POBS09.mat') && ~strcmp(POBS_dir(i).name,'POBS10.mat')
    %     continue
    % end
    
    load(['../pressure_data/' POBS_dir(i).name],'calInfoAll1')
    if isempty(calInfoAll1)
        continue
    elseif length(calInfoAll1.i0)<=4
        continue
    end
    load(['../pressure_data/' POBS_dir(i).name])

    % get indices of first point after calibration
    inorm=find(diff(data.lNormState1)==1);
    inorm=inorm+1;

    lseg=86400*7-3600; % 7 days between calibrations, minus an hour for buffer

    %----- get transients and make general plot of problem
    figure(6); clf; hold on
    for j=1:length(inorm)-1 % last segment likely too short
        % extract data segments
        transients(j).t=data.ta(inorm(j):inorm(j)+lseg);
        transients(j).p1=data.a1(inorm(j):inorm(j)+lseg);
        transients(j).p2=data.a2(inorm(j):inorm(j)+lseg);
        transients(j).T1=data.Ta1(inorm(j):inorm(j)+lseg);
        transients(j).T2=data.Ta2(inorm(j):inorm(j)+lseg);

        % correct for drift
        d1=linspace(calInfoAll1.pCal(j),calInfoAll1.pCal(j+1),length(transients(j).p1))';
        transients(j).p1_cal=transients(j).p1-d1;
        d2=linspace(calInfoAll2.pCal(j),calInfoAll2.pCal(j+1),length(transients(j).p2))';
        transients(j).p2_cal=transients(j).p2-d2;

        % decimate to 1 sample/minute for plotting
        ttemp=transients(j).t(1:60:end);
        ptemp=transients(j).p1_cal(1:60:end)-transients(j).p2_cal(1:60:end);

        yyaxis left
        plot(ttemp,ptemp,'k-','linewidth',1)
        ylabel('P (hPa)')
        datetick('x',6)
        title([staname{k} ' Calibrated, Stitched Data'])
        set(gca,'fontsize',14)
        yyaxis right
        plot([ttemp(1) ttemp(end)],[d1(1)-d2(1) d1(end)-d2(end)],'ro-','linewidth',1)
        ylabel('Calibrations (hPa)')
        box on; grid on
    end

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    % print(['../figures/stitching/stitch_assessment/transient_stacks/' staname{k} '/' staname{k}],'-dpng','-r300')

    % remove anomalous segments (empirically determined)
    if strcmp(POBS_dir(i).name,'POBS1.mat')
        transients([42 11 41 14 12])=[];
    elseif strcmp(POBS_dir(i).name,'POBS10.mat')
        transients([19 18])=[];
    elseif strcmp(POBS_dir(i).name,'POBS16.mat')
        transients(48)=[];
    elseif strcmp(POBS_dir(i).name,'POBS2.mat')
        transients([5 12 15 27])=[];
    elseif strcmp(POBS_dir(i).name,'POBS3.mat')
        transients([35 30])=[];
    elseif strcmp(POBS_dir(i).name,'POBS4.mat')
        transients([3 47 51])=[];
    end

    %----- detrend data to better isolate transient, generate stacks and mean
    figure(13); clf;
    subplot(211);hold on
    for j=1:length(transients) % last segment likely too short
        % remove slope (based on last 5 days)
        ttemp=transients(j).t-transients(j).t(1);
        m1=polyfit(ttemp(end-86400*5:end),transients(j).p1(end-86400*5:end),1);
        transients(j).p1_detrend=transients(j).p1-polyval(m1,ttemp);
        transients(j).p1_detrend=transients(j).p1_detrend-transients(j).p1_detrend(1);
        m2=polyfit(ttemp(end-86400*5:end),transients(j).p2(end-86400*5:end),1);
        transients(j).p2_detrend=transients(j).p2-polyval(m2,ttemp);
        transients(j).p2_detrend=transients(j).p2_detrend-transients(j).p2_detrend(1);

        % decimate to 1 sample/15 minutes for plotting
        ttemp=transients(j).t(1:15*60:end)-transients(j).t(1);
        ptemp=transients(j).p1_detrend(1:15*60:end)-transients(j).p2_detrend(1:15*60:end);

        % add to plot
        plot(ttemp,ptemp,'-','linewidth',1)
        ylabel('P (hPa)')
        datetick('x',7)
        title([staname{k} ' Detrended Gauge Diff. Transients'])
        set(gca,'fontsize',14)
        box on; grid on
    end

    p1_stack=cat(2,transients.p1_detrend); p1_sum=sum(p1_stack,2)/length(transients);
    p2_stack=cat(2,transients.p2_detrend); p2_sum=sum(p2_stack,2)/length(transients);

    subplot(212); hold on
    plot(transients(1).t-transients(1).t(1),p1_sum-p2_sum,'linewidth',1)
    ylabel('P (hPa)')
    datetick('x',7)
    title([staname{k} ' Mean Detrended Gauge Diff. Transient'])
    set(gca,'fontsize',14)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    % print(['../figures/stitching/stitch_assessment/transient_stacks/' staname{k} '/' staname{k} '_detrend'],'-dpng','-r300')
    % now zoom in
    xlim([0 0.2]); datetick('x',15,'keeplimits')
    subplot(211); xlim([0 0.2]); datetick('x',15,'keeplimits')
    % print(['../figures/stitching/stitch_assessment/transient_stacks/' staname{k} '/' staname{k} '_detrend_zoom'],'-dpng','-r300')

    %----- time shift to find optimal alignment between transients
    % % REMOVED BECAUSE DOES NOT WORK WELL FOR SHALLOW INSTRUMENTS
    % % which transient starts "latest"?
    % for j=1:length(transients)
    %     ptemp=transients(j).p1_detrend-transients(j).p2_detrend;
    %     pdif(j)=ptemp(8640)-ptemp(1);
    % end
    % [~,iref]=min(abs(pdif));

    % WORKS FOR ALL INSTRUMENTS
    % which transient starts "latest"?
    iloc=[];
    for j=1:length(transients)
        ptemp=transients(j).p1_detrend-transients(j).p2_detrend;
        ptemp=ptemp(1:864);
        [~,iloc(j)]=max(ptemp);
    end
    [~,iref]=min(iloc);

    pref=transients(iref).p1_detrend(1:864)-transients(iref).p2_detrend(1:864);
    pref=pref-pref(1);

    % now use misfit of first hour to align all others to reference
    figure(12); clf;
    subplot(211);hold on
    ic=[];
    for j=1:length(transients)
        tmax=300; % allow up to 5 minute shift
        for l=0:tmax
            ptemp=transients(j).p1_detrend((1:864)+l)-transients(j).p2_detrend((1:864)+l);
            ptemp=ptemp-ptemp(1);
            
            check(l+1)=std(ptemp-pref);
        end
        [~,ic(j)]=min(check);
        ic(j)=ic(j)-1; % corrects offset from indexing
        if ic(j)==tmax
            keyboard % check larger time shifts
        end

        transients(j).t_align=transients(iref).t(1:end-tmax)-transients(iref).t(1);
        transients(j).p1_align=transients(j).p1_detrend((1:end-tmax)+ic(j));
        transients(j).p1_align=transients(j).p1_align-transients(j).p1_align(1);
        transients(j).p2_align=transients(j).p2_detrend((1:end-tmax)+ic(j));
        transients(j).p2_align=transients(j).p2_align-transients(j).p2_align(1);

        % decimate to 1 sample/15 minutes for plotting
        ttemp=transients(j).t_align(1:15*60:end);
        ptemp=transients(j).p1_align(1:15*60:end)-transients(j).p2_align(1:15*60:end);

        % add to plot
        plot(ttemp,ptemp,'-','linewidth',1)
        ylabel('P (hPa)')
        datetick('x',7)
        title([staname{k} ' Aligned, Detrended Gauge Diff. Transients'])
        set(gca,'fontsize',14)
        box on; grid on
    end
    disp('time shifts:')
    disp(ic)

    p1_stack=cat(2,transients.p1_align); p1_sum=sum(p1_stack,2)/length(transients);
    p2_stack=cat(2,transients.p2_align); p2_sum=sum(p2_stack,2)/length(transients);

    subplot(212); hold on
    plot(transients(1).t_align,p1_sum-p2_sum,'linewidth',1)
    ylabel('P (hPa)')
    datetick('x',7)
    title([staname{k} ' Mean Aligned, Detrended Gauge Diff. Transient'])
    set(gca,'fontsize',14)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    % print(['../figures/stitching/stitch_assessment/transient_stacks/' staname{k} '/' staname{k} '_detrend_align'],'-dpng','-r300')
    % now zoom in
    xlim([0 0.2]); datetick('x',15,'keeplimits')
    subplot(211); xlim([0 0.2]); datetick('x',15,'keeplimits')
    % print(['../figures/stitching/stitch_assessment/transient_stacks/' staname{k} '/' staname{k} '_detrend_align_zoom'],'-dpng','-r300')

    %----- another version of this figure that shows the standard deviation
    % set end to 0 (assumes transient recovers)
    for n=1:length(transients)
        p1_stack(:,n)=p1_stack(:,n)-mean(p1_stack(end-3*86400:end,n));
        p2_stack(:,n)=p2_stack(:,n)-mean(p2_stack(end-3*86400:end,n));
    end
    p1_sum=sum(p1_stack,2)/length(transients);
    p2_sum=sum(p2_stack,2)/length(transients);

    % also get temperature stack/mean
    T1_stack=cat(2,transients.T1); T1_sum=sum(T1_stack,2)/length(transients);
    T1_sum=T1_sum-T1_sum(1);
    T2_stack=cat(2,transients.T2); T2_sum=sum(T2_stack,2)/length(transients);
    T2_sum=T2_sum-T2_sum(1);

    pstd=std(p1_stack'-p2_stack');
    figure(4); clf
    subplot(311); hold on % full duration
    hp=plot(transients(1).t_align,p1_sum-p2_sum,'k','linewidth',1);
    hs=plot(transients(1).t_align,-pstd,'r','linewidth',1);
    plot(transients(1).t_align,pstd,'r','linewidth',1)
    ylabel('P (hPa)')
    datetick('x',7)
    title([staname{k} ' Mean Detrended Gauge Diff. Transient'])
    set(gca,'fontsize',14)
    box on; grid on
    xlabel('Days')
    yyaxis right
    hc=plot(transients(1).t-transients(1).t(1),T1_sum-T2_sum,'-','linewidth',1);
    legend([hp hs hc],'mean','+/- std','Temperature','location','southeast')
    ylabel('T (C)')
    subplot(312); hold on % approx. 5 hours
    hp=plot(transients(1).t_align,p1_sum-p2_sum,'k','linewidth',1);
    hs=plot(transients(1).t_align,-pstd,'r','linewidth',1);
    plot(transients(1).t_align,pstd,'r','linewidth',1)
    ylabel('P (hPa)')
    datetick('x',7)
    set(gca,'fontsize',14)
    box on; grid on
    xlim([0 0.2]); datetick('x',15,'keeplimits')
    xlabel('Hours')
    yyaxis right
    hc=plot(transients(1).t-transients(1).t(1),T1_sum-T2_sum,'-','linewidth',1);
    legend([hp hs hc],'mean','+/- std','Temperature','location','southeast')
    ylabel('T (C)')
    subplot(313); hold on % approx. 1 hour
    hp=plot(transients(1).t_align,p1_sum-p2_sum,'k','linewidth',1);
    hs=plot(transients(1).t_align,-pstd,'r','linewidth',1);
    plot(transients(1).t_align,pstd,'r','linewidth',1)
    ylabel('P (hPa)')
    datetick('x',7)
    set(gca,'fontsize',14)
    box on; grid on
    xlim([0 0.05]); datetick('x',15,'keeplimits')
    xlabel('Hours')
    yyaxis right
    hc=plot(transients(1).t-transients(1).t(1),T1_sum-T2_sum,'-','linewidth',1);
    legend([hp hs hc],'mean','+/- std','Temperature','location','southeast')
    ylabel('T (C)')

    pause(2) % prior figure gets saved in its place if I don't build in a pause

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../figures/stitching/stitch_assessment/transient_stacks/' staname{k} '/' staname{k} '_std'],'-dpng','-r300')

    %----- correct each segment for tides, then re-create equivalent plots
    for j=1:length(transients)
        [h1,th]=oceantide(stalat(k),stalon(k),-transients(j).t(1),-transients(j).t(end),60,'osu','datenum');
        h1s=interp1(th,h1,transients(j).t,'linear','extrap');
        h1s=(h1s-mean(h1s))*100;

        % scaling
        m1=inv(h1s'*h1s)*h1s'*transients(j).p1;
        m2=inv(h1s'*h1s)*h1s'*transients(j).p2;

        % correction
        transients(j).p1_tf=transients(j).p1-m1*h1s; transients(j).p1_tf=transients(j).p1_tf-transients(j).p1_tf(1);
        transients(j).p2_tf=transients(j).p2-m2*h1s; transients(j).p2_tf=transients(j).p2_tf-transients(j).p2_tf(1);

        % remove slope (based on last 5 days)
        ttemp=transients(j).t-transients(j).t(1);
        m1=polyfit(ttemp(end-86400*5:end),transients(j).p1_tf(end-86400*5:end),1);
        transients(j).p1_tf_detrend=transients(j).p1_tf-polyval(m1,ttemp);
        transients(j).p1_tf_detrend=transients(j).p1_tf_detrend-transients(j).p1_tf_detrend(1);
        m2=polyfit(ttemp(end-86400*5:end),transients(j).p2_tf(end-86400*5:end),1);
        transients(j).p2_tf_detrend=transients(j).p2_tf-polyval(m2,ttemp);
        transients(j).p2_tf_detrend=transients(j).p2_tf_detrend-transients(j).p2_tf_detrend(1);

        % align, using lags determined above
        transients(j).p1_tf_align=transients(j).p1_tf_detrend((1:end-tmax)+ic(j));
        transients(j).p1_tf_align=transients(j).p1_tf_align-transients(j).p1_tf_align(1);
        transients(j).p2_tf_align=transients(j).p2_tf_detrend((1:end-tmax)+ic(j));
        transients(j).p2_tf_align=transients(j).p2_tf_align-transients(j).p2_tf_align(1);
    end

    p1tf_stack=cat(2,transients.p1_tf); p1tf_sum=sum(p1tf_stack,2)/length(transients);
    p1tf_sum=p1tf_sum-mean(p1tf_sum(end-3*86400:end));
    p2tf_stack=cat(2,transients.p2_tf); p2tf_sum=sum(p2tf_stack,2)/length(transients);
    p2tf_sum=p2tf_sum-mean(p2tf_sum(end-3*86400:end));

    p1tfa_stack=cat(2,transients.p1_tf_align); p1tfa_sum=sum(p1tfa_stack,2)/length(transients);
    p1tfa_sum=p1tfa_sum-mean(p1tfa_sum(end-3*86400:end));
    p2tfa_stack=cat(2,transients.p2_tf_align); p2tfa_sum=sum(p2tfa_stack,2)/length(transients);
    p2tfa_sum=p2tfa_sum-mean(p2tfa_sum(end-3*86400:end));

    T1_sum=T1_sum-mean(T1_sum(end-3*86400:end));
    T2_sum=T2_sum-mean(T2_sum(end-3*86400:end));

    figure(14); clf
    subplot(311); hold on
    plot(transients(1).t-transients(1).t(1),p1tf_sum,'linewidth',1)
    ylabel('P (hPa)')
    datetick('x',7)
    title([staname{k} ' Mean Tidally Filtered, Detrended Gauge 1 Transient'])
    set(gca,'fontsize',14)
    xlabel('Days')
    box on; grid on
    yyaxis right
    plot(transients(1).t-transients(1).t(1),T1_sum,'linewidth',1)
    ylabel('T (C)')
    subplot(312); hold on
    plot(transients(1).t-transients(1).t(1),p1tf_sum,'linewidth',1)
    ylabel('P (hPa)')
    datetick('x',7)
    set(gca,'fontsize',14)
    xlabel('Days')
    box on; grid on
    yyaxis right
    plot(transients(1).t-transients(1).t(1),T1_sum,'linewidth',1)
    ylabel('T (C)')
    xlim([0 0.2]); datetick('x',15,'keeplimits')
    subplot(313); hold on
    plot(transients(1).t-transients(1).t(1),p1tf_sum,'linewidth',1)
    ylabel('P (hPa)')
    datetick('x',7)
    set(gca,'fontsize',14)
    xlabel('Days')
    box on; grid on
    yyaxis right
    plot(transients(1).t-transients(1).t(1),T1_sum,'linewidth',1)
    ylabel('T (C)')
    xlim([0 0.05]); datetick('x',15,'keeplimits')

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../figures/stitching/stitch_assessment/transient_stacks/' staname{k} '/' staname{k} '_detrend_tfilt_G1'],'-dpng','-r300')

    figure(15); clf
    subplot(311); hold on
    plot(transients(1).t-transients(1).t(1),p2tf_sum,'linewidth',1)
    ylabel('P (hPa)')
    datetick('x',7)
    title([staname{k} ' Mean Tidally Filtered, Detrended Gauge 2 Transient'])
    set(gca,'fontsize',14)
    xlabel('Days')
    box on; grid on
    yyaxis right
    plot(transients(1).t-transients(1).t(1),T2_sum,'linewidth',1)
    ylabel('T (C)')
    subplot(312); hold on
    plot(transients(1).t-transients(1).t(1),p2tf_sum,'linewidth',1)
    ylabel('P (hPa)')
    datetick('x',7)
    set(gca,'fontsize',14)
    xlabel('Days')
    box on; grid on
    yyaxis right
    plot(transients(1).t-transients(1).t(1),T2_sum,'linewidth',1)
    ylabel('T (C)')
    xlim([0 0.2]); datetick('x',15,'keeplimits')
    subplot(313); hold on
    plot(transients(1).t-transients(1).t(1),p2tf_sum,'linewidth',1)
    ylabel('P (hPa)')
    datetick('x',7)
    set(gca,'fontsize',14)
    xlabel('Days')
    box on; grid on
    yyaxis right
    plot(transients(1).t-transients(1).t(1),T2_sum,'linewidth',1)
    ylabel('T (C)')
    xlim([0 0.05]); datetick('x',15,'keeplimits')

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../figures/stitching/stitch_assessment/transient_stacks/' staname{k} '/' staname{k} '_detrend_tfilt_G2'],'-dpng','-r300')

    % for the aligned version, let's plot G1 and G2 on same axes
    figure(16); clf
    subplot(311); hold on
    plot(transients(1).t_align,p1tfa_sum,'b','linewidth',1)
    plot(transients(1).t_align,p2tfa_sum,'r','linewidth',1)
    ylabel('P (hPa)')
    datetick('x',7)
    title([staname{k} ' Mean Tidally Filtered, Detrended Transient'])
    set(gca,'fontsize',14)
    xlabel('Days')
    box on; grid on
    yyaxis right
    plot(transients(1).t-transients(1).t(1),T1_sum,'b:','linewidth',1)
    plot(transients(1).t-transients(1).t(1),T2_sum,'r:','linewidth',1)
    ylabel('T (C)')
    legend('p1','p2','T1','T2','location','northeast')
    subplot(312); hold on
    plot(transients(1).t_align,p1tfa_sum,'b','linewidth',1)
    plot(transients(1).t_align,p2tfa_sum,'r','linewidth',1)
    ylabel('P (hPa)')
    datetick('x',7)
    set(gca,'fontsize',14)
    xlabel('Days')
    box on; grid on
    yyaxis right
    plot(transients(1).t-transients(1).t(1),T1_sum,'b:','linewidth',1)
    plot(transients(1).t-transients(1).t(1),T2_sum,'r:','linewidth',1)
    ylabel('T (C)')
    legend('p1','p2','T1','T2','location','northeast')
    xlim([0 0.5]); datetick('x',15,'keeplimits')
    subplot(313); hold on
    plot(transients(1).t_align,p1tfa_sum,'b','linewidth',1)
    plot(transients(1).t_align,p2tfa_sum,'r','linewidth',1)
    ylabel('P (hPa)')
    datetick('x',7)
    set(gca,'fontsize',14)
    xlabel('Days')
    box on; grid on
    yyaxis right
    plot(transients(1).t-transients(1).t(1),T1_sum,'b:','linewidth',2)
    plot(transients(1).t-transients(1).t(1),T2_sum,'r:','linewidth',2)
    ylabel('T (C)')
    legend('p1','p2','T1','T2','location','northeast')
    xlim([0 0.01]); datetick('x',15,'keeplimits')

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../figures/stitching/stitch_assessment/transient_stacks/' staname{k} '/' staname{k} '_detrend_align_tfilt_G1G2'],'-dpng','-r300')

    % % unless I decimate the data, these files are just too large to save
    % save(['../transients/' staname{k}],'transients')
    clearvars transients
end
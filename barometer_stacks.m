% barometer_stacks.m
%
% Stacked plots of barometer curves and ideal gas law equivalents
%

clear; close all

figure(3); clf; hold on
figure(4); clf
n=0;

POBS_dir=dir('../pressure_data/');
POBS_list={POBS_dir.name}';
file_check=cellfun(@(v)v(1),POBS_list);
i_list=find(~eq(file_check,'.'));
for i=1:length(i_list)
    ii=i_list(i);

    % load barometer and calibration files
    load([POBS_dir(ii).folder '/' POBS_dir(ii).name],'barInfoAll','calInfoAll1','calInfoAll2')
    if isempty(barInfoAll)
        continue
    elseif length(barInfoAll.i0p)<=4
        continue
    end

    % convert filename to label
    if length(POBS_dir(ii).name)==9
        snum=str2double(POBS_dir(ii).name(5));
    elseif length(POBS_dir(ii).name)==10
        snum=str2double(POBS_dir(ii).name(5:6));
    else
        snum=str2double(POBS_dir(ii).name(5));
    end
    sname=['POBS-' num2str(snum)];

    % start plotting
    n=n+1;

    figure(3)
    t=barInfoAll.t0p;
    if ~strcmp(POBS_dir(ii).name,'POBS2.mat')
        p1=barInfoAll.pCal;
        if strcmp(POBS_dir(ii).name,'POBS1.mat')
            t=t(1:end-5); p1=p1(1:end-5);
        elseif strcmp(POBS_dir(ii).name,'POBS15.mat')
            t=t(1:end-10); p1=p1(1:end-10);
        end
        h1=plot(t,p1-mean(p1)+(n-1)*10,'ob','markersize',10,'linewidth',1);
    end
    p2=barInfoAll.T;
    if strcmp(POBS_dir(ii).name,'POBS1.mat')
        p2=p2(1:end-5);
    elseif strcmp(POBS_dir(ii).name,'POBS15.mat')
        p2=p2(1:end-10);
    end
    p2=(p2-mean(p2))*3.71;
    h2=plot(t,p2+(n-1)*10,'^r','markersize',10,'linewidth',1);

    if strcmp(POBS_dir(ii).name,'POBS09.mat')
        text(t(1)-75,mean(p2)+(n-1)*10,sname,'fontsize',12)
    else
        text(t(end)+15,mean(p2)+(n-1)*10,sname,'fontsize',12)
    end

    figure(4)
    t=barInfoAll.t0p;
    if ~strcmp(POBS_dir(ii).name,'POBS2.mat')
        p1=barInfoAll.pCal;
        if strcmp(POBS_dir(ii).name,'POBS1.mat')
            t=t(1:end-5); p1=p1(1:end-5);
        elseif strcmp(POBS_dir(ii).name,'POBS15.mat')
            t=t(1:end-10); p1=p1(1:end-10);
        end
        h1=plot(t,p1-mean(p1)+(n-1)*10,'ob','markersize',10,'linewidth',1);
    end
    p2=barInfoAll.T;
    if strcmp(POBS_dir(ii).name,'POBS1.mat')
        p2=p2(1:end-5);
    elseif strcmp(POBS_dir(ii).name,'POBS15.mat')
        p2=p2(1:end-10);
    end
    p2=(p2-mean(p2))*3.71;
    h2=plot(t,p2+(n-1)*10,'^r','markersize',10,'linewidth',1);

    if strcmp(POBS_dir(ii).name,'POBS09.mat')
        text(t(1)-75,mean(p2)+(n-1)*10,sname,'fontsize',12)
    else
        text(t(end)+15,mean(p2)+(n-1)*10,sname,'fontsize',12)
    end
end

figure(3)
ylim([-10 130])
ylabel('P (hPa)')
legend('Barometer','Ideal Gas Law')
datetick('x',6)
set(gca,'fontsize',14)
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures/manuscript/supplement/barometer_stack','-dpng','-r300')
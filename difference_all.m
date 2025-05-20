% difference_all.m

clear; close all

geom=load('../pressure_data/geometry.mat');

POBS_dir=dir('../stitched_data/drift_corrected');
POBS_list={POBS_dir.name}';
file_check=cellfun(@(v)v(1),POBS_list);
i_list=find(eq(file_check,'P')); % subselects only files begining with 'P'

t_hf=[];
p_hf=[];
staname=[];
for i=1:length(i_list)
    k=i_list(i);
    load([POBS_dir(k).folder '/' POBS_dir(k).name],'dataf')
    t_hf{i}=dataf.tf;
    p_hf{i}=dataf.p2_dcor;
    d(i)=round(mean(dataf.p2_dcor)/100);
    staname{i}=POBS_dir(k).name(1:end-4);
    if length(staname{i})>6
        staname{i}=staname{i}(1:5);
    end
    j=find(strcmp(staname{i},geom.filename));
    if strcmp(staname{i},'POBS7')
        j=6;
    end
    stalon(i)=geom.stalon(j);
    stalat(i)=geom.stalat(j);
    stadepth(i)=geom.stadepth(j);
end
n=length(t_hf);

G1=load('../../hikurangi/processed_data/GONDOR-I.mat');
for i=1:length(G1.pf)
    t_hf{n+i}=G1.tf{i}';
    p_hf{n+i}=G1.pf{i}';
    d(n+i)=G1.stadepth(i);
    staname{n+i}=G1.staname{i};
    stalon(n+i)=G1.stalon(i);
    stalat(n+i)=G1.stalat(i);
    stadepth(n+i)=G1.stadepth(i);
end

%% network map (A-0-A and APG)
addpath('../../hikurangi/code/m_map')

% setup base map
figure(54); clf; hold on
lat1=-40; lat2=-38;
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
m_plot(stalon(1:n),stalat(1:n),'^k','linewidth',1,'markerfacecolor','y','markersize',10)
m_text(stalon(1:n)-0.25,stalat(1:n)-0.02,staname(1:n),'fontsize',12)

m_plot(stalon(n+1:end),stalat(n+1:end),'sk','linewidth',1,'markerfacecolor','c','markersize',10)
m_text(stalon(n+1:end)-0.25,stalat(n+1:end)-0.02,staname(n+1:end),'fontsize',12)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures/differences/combined_map','-dpng','-r300')
print('../figures/differences/combined_map','-depsc','-vector')

%% differences

std_mat=zeros(length(t_hf));
for j=1:length(t_hf)
    t1=t_hf{j};
    p1=p_hf{j}-mean(p_hf{j});
    
    for k=1:length(t_hf)
        if k==j
            continue
        end
        t2=t_hf{k};
        p2=p_hf{k}-mean(p_hf{k});

        [tdif,ia,ib]=intersect(round(t1,6),round(t2,6));
        pdif=p1(ia)-p2(ib);
        if j>n || k>n % detrend non-A0A differences
            pdif=detrend(pdif);
        end

        std_mat(j,k)=std(detrend(pdif));

        % figure(2); clf; hold on
        % plot(tdif,p1(ia),'linewidth',0.5)
        % plot(tdif,p2(ib),'linewidth',0.5)
        % plot(tdif,pdif,'k','linewidth',1)
        % legend(staname{j},staname{k},'diff','location','best')
        % ylabel('\DeltaP (cm)')
        % datetick('x',6)
        % title([staname{j} ' (' num2str(d(j)) ' m) -- ' staname{k} ' (' num2str(d(k)) ' m)'])
        % set(gca,'fontsize',14)
        % box on; grid on
        % 
        % fh=gcf;
        % fh.PaperUnits='inches';
        % fh.PaperPosition=[0 0 11 8.5];
        % print(['../figures/differences/all/' staname{j} '-' staname{k}],'-djpeg','-r100')
    end
end

std_mat(isnan(std_mat))=0;
srt_std=sort(std_mat(:));

imax=length(srt_std):-2:length(srt_std)-100;

for i=1:length(imax)
    [a,b]=find(std_mat==srt_std(imax(i)));
    disp(['Stations ' staname{a(1)} ' (' num2str(stadepth(a(1))) ') and ' ...
        staname{b(1)} ' (' num2str(stadepth(b(1))) ') poorly correlated'])
end
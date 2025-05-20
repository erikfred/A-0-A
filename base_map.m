% base_map.m
%
% coastline, bathymetry, etc. upon which to plot various deployment
% geometries
%

clear; close all

addpath('../../hikurangi/code/m_map')

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

%% change this section to mark different things on the map

% overlapping GONDOR-I and -II A-0-A locations
g1=load('../pressure_data/geometry.mat');
g2=load('../pressure_data_Y2/geometry.mat');

% Y1 markers
m_plot(g1.stalon,g1.stalat,'^k','linewidth',1,'markerfacecolor','c','markersize',10)
m_text(g1.stalon-0.21,g1.stalat,g1.staname,'fontsize',12)

% Y2 markers
m_plot(g2.stalon,g2.stalat,'sk','linewidth',1,'markerfacecolor','y','markersize',10)
m_text(g2.stalon+0.1,g2.stalat,g2.staname,'fontsize',12,'color','c')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures_Y2/maptest','-dpng','-r300')
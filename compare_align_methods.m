% compare_align_methods.m
%
% Make stem plots to demonstrate how slopes change with various alignment
% methods. Make separate plots to show how slopes improve yet more after
% correction with ocean models.
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

load('../figures_Y2/combine_years/comb_rates.mat')
load('../figures_Y2/combine_years/Y1_rates.mat')

% do some pruning
a2r(1)=a1r(1); b2r(1)=b1r(1); % POBS-04 has bad data on gauge 2
p9rate2(1)=p9rate1(1); % POBS-04 has bad data on gauge 2
cpierate2(1)=cpierate1(1); % POBS-04 has bad data on gauge 2
glorate2(1)=glorate1(1); % POBS-04 has bad data on gauge 2
gcorrate2(1)=gcorrate1(1); % POBS-04 has bad data on gauge 2

sname([5 12])=[]; % POBS-10 and POBS-09 aren't relevant here
a2r([5 12])=[]; b2r([5 12])=[]; % POBS-10 and POBS-09 aren't relevant here

% slopes before and after drift correction
figure(7); clf; hold on
h1=stem((1:10)+0.1,a2r,'ko','markersize',12,'markerfacecolor','r');
h2=stem((1:10)+0.1,b2r','ko','markersize',12,'markerfacecolor',[0 114 189]/255);
% view(90,-90)
ylabel('Slope (hPa/yr)')
legend('Original','Drift-Corrected','location','southeast')
set(gca,'fontsize',20)
set(gca,'xtick',1:12)
set(gca,'xticklabels',sname)
xlim([0.1 10.9])
box on; grid on
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/presentations/trend_summary1','-dpng','-r300')
print('../figures/presentations/trend_summary1','-depsc','-vector')

% manual entry of Y1 model-corrected rates
y1r=[3.7 4.2 3.5 5.1 0.6 2.2 0.0 4.0 0.9 1.9];

delete(h1)
h11=stem((1:10)+0.1,y1r,'vk','markersize',12,'markerfacecolor','w');
legend('Drift-Corrected','Drift+Model-Corrected','location','southeast')
ylim([-5 10])
print('../figures/presentations/trend_summary2','-dpng','-r300')
print('../figures/presentations/trend_summary2','-depsc','-vector')

delete(h11)
h3=stem((1:10)+0.1,p9rate2,'sk','markersize',12,'markerfacecolor','g');
h4=stem((1:10)+0.1,cpierate2,'^k','markersize',12,'markerfacecolor','m');
h5=stem((1:10)+0.1,glorate2,'<k','markersize',12,'markerfacecolor','y');
legend('Year 1 Drift-Corrected','POBS-09 Align','CPIES Align','GLORYS Align','location','southeast')
ylim([-5 10])
print('../figures/presentations/trend_summary3','-dpng','-r300')
print('../figures/presentations/trend_summary3','-depsc','-vector')

delete(h2); delete(h3); delete(h4); delete(h5)
h6=stem((1:10)+0.1,y1r,'vk','markersize',12,'markerfacecolor','w');
h7=stem((1:10)+0.1,gcorrate2,'pk','markersize',15,'markerfacecolor','k');
legend('Year1 Model-Corrected','Both Years Model-Corrected','location','southeast')
ylim([-4 6])
print('../figures/presentations/trend_summary4','-dpng','-r300')
print('../figures/presentations/trend_summary4','-depsc','-vector')

% predicted rates from TDEFNODE
tr2=[0.4 0.7 0.8 0.0 0.4 0.0 0.5 0.0 0.1 0.6];

delete(h6)
h8=stem((1:10)+0.1,tr2,'xk','markersize',15,'linewidth',1);
legend('Both Years Model-Corrected','Inter-SSE Model','location','southeast')
ylim([-3 3])
print('../figures/presentations/trend_summary5','-dpng','-r300')
print('../figures/presentations/trend_summary5','-depsc','-vector')

delete(h7)
legend('Inter-SSE Model','location','southeast')
ylim([-3 3])
print('../figures/presentations/trend_summary0','-dpng','-r300')
print('../figures/presentations/trend_summary0','-depsc','-vector')
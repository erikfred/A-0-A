% geometry_Y2.m
%
% read file/geometry information and save as .mat for easy use
%

clear; close all

y2dir=dir('../../hikurangi/processed_POBS_Y2');
flist={y2dir.name}';
fcheck=cellfun(@(v)v(1),flist);
i_list=find(eq(fcheck,'P'));
filename=flist(i_list)';

for i=1:length(filename)
    staname{i}=[filename{i}(1:4) '-' filename{i}(5:end)];
end

D=readtable('../POBS_locations_Y2.csv');
stalat=table2array(D(:,3));
stalon=table2array(D(:,4));
stadepth=table2array(D(:,5));

save('../pressure_data_Y2/geometry','filename','staname','stalat','stalon','stadepth')
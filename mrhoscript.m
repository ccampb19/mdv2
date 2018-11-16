% mrhoscript
clear
close all

opt.basedir = 'D:\fdpm.data\tempdrift7\181030\';
opt.nind = 1.332; %1.33 for water, 1.36 for 20% intralipid, 1.4 for PDMS  
opt.laser_names=[660,688,781,806,828,849];  %only for plotting names
opt.usediodes = 1:6;

opt.geomadjust = 1;

%Load all possible rhos, to be adjusted in endmat
opt.rhorange = 11:30;

% Enter the first filename to be loaded, less the '-dcswitch.asc'
% Code will replace the rho in the prototype with opt.rhorange.
% Ex: 'ex-12-1-baseline' will load 'ex-12-1-baseline-dcswitch.asc' thru 
%     'ex-30-1-baseline-dcswitch.asc'
opt.filenameprototype = '11-1-mdfinal';
opt.darkname = 'dark';
opt.darkreps = 3;

%Startrhos must be at least min(rhorange),
%and probably don't need to change much.
%Endrhos must be at least max(rhorange)-2.
%Endfreqs don't need to be exact. 

endmat = [
% [startrho,endrho,endfreq]
    12 28 600;         % 660
    12 28 600;         % 690
    12 28 600;         % 785
    12 28 600;         % 810
    12 28 600;         % 830
    12 28 600;         % 850
    ];

%% Do the fit
data = mdload(opt);
outdata = mdprocess(opt,endmat,data);
mdplot(opt,outdata);
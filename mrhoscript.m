% mrhoscript
clear
% close all

opt.basedir = 'D:\fdpm.data\troubleshooting\190130\';
opt.nind = 1.332; %1.33 for water, 1.36 for 20% intralipid, 1.4 for PDMS  
opt.laser_names=[660,688,781,806,828,849];  %only for plotting names
opt.usediodes = 1:6;

% Adjust rho for LBS fiber positioning in ferrules?
opt.geomadjust = 0;

% Fit Broadband? One sphere file only (full filename)
opt.bb = 1;
opt.sphname = 'sphere-tis.asc';

% Load all possible rhos, to be adjusted in endmat
startrho = 11;
endrho = 28;
opt.rhosteps = 1;
opt.rhorange = startrho:opt.rhosteps:endrho;

% Enter the first filename to be loaded, less the '-dcswitch.asc'
% Code will replace the rho in the prototype with opt.rhorange.
% Ex: 'ex-12-1-baseline' will load 'ex-12-1-baseline-dcswitch.asc' thru 
%     'ex-30-1-baseline-dcswitch.asc'
opt.filenameprototype = [num2str(startrho) '-1il_only_10460_1000_noled'];

% fdpm dark msmts
opt.darkname = 'dark';
opt.darkreps = 5;

% Startrhos must be at least min(rhorange),
% and probably don't need to change much.
% Endrhos must be at least max(rhorange)-2.
% Endfreqs don't need to be exact. 

%% Do the fit
endmat = repmat([startrho,endrho,560],6,1);

data = mdload(opt);
tic
outdata = mdprocess(opt,endmat,data);
toc
phantom = mdplot(opt,outdata);
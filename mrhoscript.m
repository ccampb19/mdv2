% mrhoscript
% clear
% close all

opt.basedir = 'D:\fdpm.data\MD-2019\190521\';
opt.nind = 1.4; %1.33 for water, 1.36 for 20% intralipid, 1.4 for PDMS  
opt.laser_names=[660,688,781,806,828,849];  %only for plotting names
opt.usediodes = 1:6;

% Adjust rho for LBS fiber positioning in ferrules?
opt.geomadjust = 1;

% Fit Broadband? One sphere file only (full filename)
opt.bb = 1;
opt.sphname = 'sphere-tis.asc';
opt.smooth = 1;
opt.threshold = 4000; % Counts where uncalibrated spectrometer loses linearity

% Load all possible rhos, to be adjusted in endmat (use long rhos if
% possible)
startrho = 12;
% With more dynamic range, broadband needs lower rhos for adequate counts
bbstartrho = 12;
endrho = 30;
opt.rhosteps = 1;
opt.rhorange = startrho:opt.rhosteps:endrho;
opt.bbrhorange = bbstartrho:opt.rhosteps:endrho;
endfreq = 560;

% Enter the first filename to be loaded, less the '-dcswitch.asc'
% Code will replace the rho in the prototype with opt.rhorange.
% Ex: 'ex-12-1-baseline' will load 'ex-12-1-baseline-dcswitch.asc' thru 
%     'ex-30-1-baseline-dcswitch.asc'
base = 'ph017';
opt.filenameprototype = [num2str(startrho) '-1-' base ];

% fdpm dark msmts
opt.darkname = 'dark';
opt.darkreps = 6;

% Startrhos must be at least min(rhorange),
% and probably don't need to change much.
% Endrhos must be at least max(rhorange)-2.
% Endfreqs don't need to be exact. 

%% Do the fit
endmat = repmat([startrho,endrho,endfreq],6,1);

data = mdload(opt);
tic
outdata = mdprocess_fd(opt,endmat,data);
if opt.bb==1
    outdata = mdprocess_bb(opt,data,outdata);
end
toc
[phantom,bbphantom] = mdplot(opt,outdata);
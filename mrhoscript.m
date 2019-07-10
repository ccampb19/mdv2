% mrhoscript
clear
% close all

opt.basedir = 'D:\fdpm.data\bovine_deox\190709\';
opt.nind = 1.33; %1.33 for water, 1.36 for 20% intralipid, 1.4 for PDMS  
opt.laser_names=[630,660,688,781,828,848];  %only for plotting names
opt.usediodes = 1:6;

% Adjust rho for LBS fiber positioning in ferrules?
% 2/27/19 Rotated ferrule
opt.geomadjust = 1;

% Fit Broadband? One sphere file only (full filename)
opt.bb = 1;
opt.bbdark = 1;
opt.sphname = 'sphere';
opt.sphreps = -1;
opt.smooth = 1;
opt.threshold = 6000; % Counts where uncalibrated spectrometer loses linearity

% Load all possible rhos, to be adjusted in endmat (use long rhos if
% possible)
startrho = 12;
% With more dynamic range, broadband needs lower rhos for adequate counts
bbstartrho = 11;
endrho = 28;
opt.rhosteps = 1;
opt.rhorange = startrho:opt.rhosteps:endrho;
opt.bbrhorange = bbstartrho:opt.rhosteps:endrho;
endfreq = 560;

% Enter the first filename to be loaded, less the '-dcswitch.asc'
% Code will replace the rho in the prototype with opt.rhorange.
% Ex: 'ex-12-1-baseline' will load 'ex-12-1-baseline-dcswitch.asc' thru 
%     'ex-30-1-baseline-dcswitch.asc'
base = 'blood_hio2_good';
opt.filenameprototype = [num2str(startrho) '-1-' base ];

% fdpm dark msmts
opt.darkname = 'dark';
opt.darkreps = 5;

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
[phantom,bbphantom,albphantom] = mdplot(opt,outdata);
% mrhoscript
clear
% close all

opt.basedir = '~/fdpm.data/hemozoin/190923';
opt.nind = 1.33; %1.33 for water, 1.36 for 20% intralipid, 1.4 for PDMS  
opt.laser_names=[631 660 689 782 828 849];  %only for plotting names
opt.usediodes = 1:6;

% Adjust rho for LBS fiber positioning in ferrules?
% 2/27/19 Rotated ferrule
opt.geomadjust = 1;

% Fit Broadband? One sphere file only (full filename)
opt.bb = 1;
opt.subtractdark = 1;
opt.overx = 1;
opt.sphname = 'sphere';
opt.sphreps = -1;
opt.smooth = 0;
opt.threshold = 6000; % Counts where uncalibrated spectrometer loses linearity
opt.binning = 6;

% Load all possible rhos, to be adjusted in endmat (use long rhos if
% possible)
startrho = 14;
% With more dynamic range, broadband needs lower rhos for adequate counts
% rhocal = linspace(11,29,20);
bbstartrho = 14;
endrho = 28;
opt.rhosteps = 1;
opt.rhorange = startrho:opt.rhosteps:endrho;
opt.bbrhorange = bbstartrho:opt.rhosteps:endrho;
endfreq = 560;

% Enter the first filename to be loaded, less the '-dcswitch.asc'
% Code will replace the rho in the prototype with opt.rhorange.
% Ex: 'ex-12-1-baseline' will load 'ex-12-1-baseline-dcswitch.asc' thru 
%     'ex-30-1-baseline-dcswitch.asc'
% base = 'blood2_DO34_reset';
base = 'blood2_DO34_reset';
opt.filenameprototype = [num2str(startrho) '-' base ];
% opt.filenameprototype = [base '-' num2str(startrho)];

% fdpm dark msmts
opt.darkname = 'dark2';
opt.darkreps = 5;

% Startrhos must be at least min(rhorange),
% and probably don't need to change much.
% Endrhos must be at least max(rhorange)-2.
% Endfreqs don't need to be exact. 

%% Do the fit
endmat = repmat([startrho,endrho,endfreq],length(opt.laser_names),1);

data = mdload(opt);
tic
outdata = mdprocess_fd(opt,endmat,data);
if opt.bb==1
    outdata = mdprocess_bb(opt,data,outdata);
end
toc
[phantom,bbphantom,albphantom] = mdplot(opt,outdata);

% Implement chromophore fit
uchrom = [1 1 1 0 1 1 1];
[cnames,fitt,concs,wvout] = mdchromfit(albphantom,'chromophores_moo_hz.txt',uchrom);
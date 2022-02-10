% mrhoscript
clear
% close all

opt.basedir = 'D:\fdpm.data\training_sn\220204';
opt.systemname = 'dcswitch'; % 'miniLBS', or 'dcswitch', for data file names
opt.nind = 1.4; %1.33 for water, 1.36 for 20% intralipid, 1.4 for PDMS  
opt.laser_names=[660 689 782 808 828 849];  %only for plotting names
opt.usediodes = 1:numel(opt.laser_names);

% Adjust rho for LBS fiber positioning in ferrules?
% 2/27/19 Rotated ferrule
opt.geomadjust = 1;

% Fit Broadband? One sphere file only (full filename)
opt.bb = 1;
opt.subtractbbdark = 1; % Dark measurements should be required, but go nuts.
opt.overx = 1; % If an overexposed tis file was created, stitch them together (optionally)
opt.sphname = 'rc04-1';
opt.sphreps = -1; % Sphere cal not recommended. -1 to skip.
opt.smooth = 0; % Smooth broadband reflectance?
opt.threshold = 7500; % Counts where uncalibrated spectrometer apparently loses linearity
opt.binning = 4;    % ~2 pixels per nm, but 9nm resolution with this spec. Binning helps smooth data

% Load all possible rhos, to be adjusted in endmat (use long rhos if
% possible)
startrho = 14;
% With more dynamic range, broadband needs lower rhos for adequate counts
bbstartrho = 14;
endrho = 28;
opt.rhosteps = 1;
opt.rhorange = startrho:opt.rhosteps:endrho;
opt.bbrhorange = bbstartrho:opt.rhosteps:endrho;

% NEW: Software writes the same filename regardless, so if broadband fiber
% is offset, set it here (-1.0 means 1mm shorter SDS than FD SDS)
opt.bboffset = -1.1;

% opt.rhorange = [11,12,17,18];
% opt.bbrhorange = opt.rhorange;
endfreq = 420;

% Enter the first filename to be loaded, except the '-dcswitch.asc'
% Code will replace the rho in the prototype with opt.rhorange.
% Ex: 'ex-12-1-baseline' will load 'ex-12-1-baseline-dcswitch.asc' thru 
%     'ex-30-1-baseline-dcswitch.asc'
base = 'ph017';
opt.filenameprototype = [num2str(startrho), '-', base];

% fdpm dark msmts
% comment out if none available
opt.darkname = 'dark';
opt.darkreps = 5;

% Startrhos must be at least min(rhorange),
% and probably don't need to change much.
% Endrhos must be at least max(rhorange)-2.
% Endfreqs don't need to be exact. 

%% Do the fit
endmat = repmat([startrho,endrho,endfreq],numel(opt.laser_names),1);

data = mdload(opt);
tic
outdata = mdprocess_fd(opt,endmat,data);
if opt.bb==1
    outdata = mdprocess_bb(opt,data,outdata);
end
toc
[phantom,bbphantom,albphantom] = mdplot(opt,outdata);

%% Chromophore fit
% Our only lipid chromophore is from rendered pig fat, I think, which
% doesn't match the phospholipids in Intralipid, so I don't use that one.
% uchrom = [1 1 1 0 1 1];
% [fitt,concs,wvout] = mdchromfit(albphantom,'chromophores_moo.txt',uchrom);
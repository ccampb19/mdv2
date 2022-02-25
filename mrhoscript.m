% mrhoscript
% be sure to have geomadjust.m file open/run in MATLAB prior to
clear
close all

opt.basedir = 'G:\My Drive\SiPM_Project\fdpm.data\220219-LBS\220219';
opt.systemname = 'dcswitch'; % 'miniLBS', or 'dcswitch', for data file names
opt.nind = 1.4; %1.33 for water, 1.36 for 20% intralipid, 1.4 for PDMS  
opt.laser_names=[660 689 782 808 828 849];  %only for plotting names reference dcswitch file
opt.usediodes = 1:numel(opt.laser_names);

% Adjust rho for LBS fiber positioning in ferrules?
% 2/27/19 Rotated ferrule
opt.geomadjust = 1;

% Fit Broadband? One sphere file only (full filename)
opt.bb = 1; % must use to get broadband ops
opt.subtractbbdark = 1; % Dark measurements should be required, but go nuts.
% overx will subtract out the dark. Normalizes everything according to 
% integration time. Turn off if weird in broadband, set to 0 to turn off
opt.overx = 0; % If an overexposed tis file was created, stitch them together (optionally)
opt.sphname = 'rc04-1';
opt.sphreps = -1; % Sphere cal not recommended. -1 to skip.
opt.smooth = 0; % Smooth broadband reflectance using moving average smoothing
opt.threshold = 7500; % Counts where uncalibrated spectrometer apparently loses linearity
opt.binning = 4;    % ~2 pixels per nm, but 9nm resolution with this spec. Binning helps smooth data

% Load all possible rhos, to be adjusted in endmat (use long rhos if
% possible)
startrho = 14; %start SDS processing at 14mm 
% With more dynamic range, broadband needs lower rhos for adequate counts
bbstartrho = 14; 
endrho = 28;
opt.rhosteps = 1;
opt.rhorange = startrho:opt.rhosteps:endrho;
opt.bbrhorange = bbstartrho:opt.rhosteps:endrho;

% NEW: Software writes the same filename regardless, so if broadband fiber
% is offset, set it here (-1.0 means 1mm shorter SDS than FD SDS)
% if the ferrules are switched and the probe is stationary then set this to
% 0 and the 
opt.bboffset = -1.1; %used because the fibers for laser and bb are 1.1mm apart

% opt.rhorange = [11,12,17,18];
% opt.bbrhorange = opt.rhorange;
endfreq = 420;

% Enter the first filename to be loaded, except the '-dcswitch.asc'
% Code will replace the rho in the prototype with opt.rhorange.
% Ex: 'ex-12-1-baseline' will load 'ex-12-1-baseline-dcswitch.asc' thru 
%     'ex-30-1-baseline-dcswitch.asc'
base = 'NorthwellB'; %case sensitive
opt.filenameprototype = [num2str(startrho), '-', base];

% fdpm dark msmts
% comment out if none available
opt.darkname = 'dark';
opt.darkreps = 10; %input number of dark measurements taken

% Startrhos must be at least min(rhorange),
% and probably don't need to change much.
% Endrhos must be at least max(rhorange)-2.
% Endfreqs don't need to be exact. 

%% Do the fit
endmat = repmat([startrho,endrho,endfreq],numel(opt.laser_names),1); %creates a matrix

data = mdload(opt);
tic
outdata = mdprocess_fd(opt,endmat,data);
if opt.bb==1
    outdata = mdprocess_bb(opt,data,outdata);
end
toc
[phantom,bbphantom,albphantom] = mdplot(opt,outdata); %plots results

abs = outdata.opfd(1,:); 
scat = outdata.opfd(2,:);

% opTable = 0;
% for i = 1:length(opt.laser_names)
%     opTable(i,1) = opt.laser_names(i);
%     opTable(i,2) = abs(i);
%     opTable(i,3) = scat(i);
% end

% this does the same as the for loop above
opTable = [ opt.laser_names', outdata.opfd(1,:)', outdata.opfd(2,:)' ];

format default %formats table so that it prints without shifting decimal points
opTable %print out the OPs in the command window, format-> 'Wavelength Abs Scat'

% write information to file 
dlmwrite('northwellB_bb.txt', bbphantom, '\t' );
dlmwrite('northwellB.txt', opTable, '\t' );


%% Chromophore fit
% Our only lipid chromophore is from rendered pig fat, I think, which
% doesn't match the phospholipids in Intralipid, so I don't use that one.
% uchrom = [1 1 1 0 1 1];
% [fitt,concs,wvout] = mdchromfit(albphantom,'chromophores_moo.txt',uchrom);
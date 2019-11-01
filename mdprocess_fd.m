function [OUTDATA] = mdprocess_fd(OPTS,ENDMAT,DATAS)
%MDLOAD loads files for multirho fits
%   OPTS contains all necessary data to load a mutlidistance dataset.
%   ENDMAT contains other information (fix descriptions)
%   DATAS contains that dataset.
%
%   PROCESSED is processed data

% Ensure ENDMAT is the correct size and in bounds
assert((size(ENDMAT,1) == size(DATAS.complex,3)) & (size(ENDMAT,2) == 3),...
    'Endmat must be size (num diodes) x 3');

endmatlength = size(ENDMAT,1);
assert(sum(ENDMAT(:,1) >= OPTS.rhorange(1))== endmatlength...
    & sum(ENDMAT(:,2) <= OPTS.rhorange(end))== endmatlength...
    & sum(ENDMAT(:,2) > ENDMAT(:,1)+2)== endmatlength...
    & sum(ENDMAT(:,3) <= DATAS.freqs(end))== endmatlength,...
    'Endmat parameters out of range');

% Correct endmat and rhos, if necessary
if isfield(OPTS,'darkname')
%     highcut = DATAS.cutoff(:,1) > OPTS.rhorange(end-2);
    ENDMAT(:,2:3) = DATAS.cutoff;
%     ENDMAT(highcut,2) = OPTS.rhorange(end-2);
end

if OPTS.geomadjust
    newrhos = geomadjust(OPTS.usediodes,OPTS.rhorange);
else
    newrhos = repmat(OPTS.rhorange,length(OPTS.usediodes),1);
end

% Do the fit one diode at a time.
% lsqoptions = optimoptions('lsqcurvefit','Algorithm','Levenberg-Marquardt',...
%     'FunctionTolerance',5e-12,'StepTolerance',1e-7,...
%     'MaxFunctionEvaluations',60000,'MaxIterations',60000);
psoptions = optimoptions('patternsearch','FunctionTolerance',1e-14,...
    'MaxFunctionEvaluations',60000,'MaxIterations',60000,...
    'MeshTolerance',2*eps);

disp('Exit Flags:')

% Where the magic happens
OUTDATA = procdiodes(ENDMAT,DATAS,OPTS,newrhos,psoptions);

OUTDATA.fderr = squeeze(std(OUTDATA.rmu,[],1));

% Fit scattering parameters
x = OPTS.laser_names(OPTS.usediodes);
pwrlaw = @(b,x) b(1).*(x.^-b(2));
fdmusp = OUTDATA.opfd(2,:);

options = optimset('MaxFunEvals',10000,'Display','off');
nrmrsd = @(b) norm(fdmusp(~isnan(fdmusp)) - pwrlaw(b,x(~isnan(fdmusp))));
OUTDATA.pwrfit = fminsearch(nrmrsd,[4000,1.3],options);
% OUTDATA.pwrfit = fminsearch(nrmrsd,[1,10,1000],options);
function [OUTDATA] = mdprocess(OPTS,ENDMAT,DATAS)
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
assert(sum(ENDMAT(:,1) >= OPTS.rhorange(1))==endmatlength...
    & sum(ENDMAT(:,2) < OPTS.rhorange(end)-1)==endmatlength...
    & sum(ENDMAT(:,2) > ENDMAT(:,1)+2)==endmatlength...
    & sum(ENDMAT(:,3) <= DATAS.freqs(end))==endmatlength,...
    'Endmat parameters out of range');

% Correct endmat and rhos, if necessary
if isfield(OPTS,'darkname')
%     highcut = DATAS.cutoff(:,1) > OPTS.rhorange(end-2);
    ENDMAT(:,2:3) = DATAS.cutoff;
%     ENDMAT(highcut,2) = OPTS.rhorange(end-2);
end

if OPTS.geomadjust
    newrhos = geomadjust(OPTS.usediodes,OPTS.rhorange);
end

% Do the fit one diode at a time.
options = optimoptions('lsqcurvefit','Algorithm','Levenberg-Marquardt',...
    'FunctionTolerance',1e-9,'StepTolerance',1e-9,...
    'MaxFunctionEvaluations',25000,'MaxIterations',25000);

endfreqidxs = zeros(6,1);
startrhoidxs = endfreqidxs;
for didx = 1:length(OPTS.usediodes)
    rhorange = newrhos(didx,:);
    fititer = 0;
    for endidx = (ENDMAT(didx,2)-5:ENDMAT(didx,2))-ENDMAT(didx,1)+1

        % Round estimated endfreqs to existing frequencies
        endfreqidxs(didx) = find(DATAS.freqs<ENDMAT(didx,3),1,'last');
        startrhoidxs(didx) = find(OPTS.rhorange == ENDMAT(didx,1));
        % For reshape operation to work as I would like, must do (elementwise)
        % transpose here.
        diodedata = DATAS.complex(:,1:endfreqidxs(didx),didx).';
        YDATA = [];
        for startidx = (startrhoidxs(didx)+1):(endidx-1)
    %         for endidx = startidx+1:(ENDMAT(didx,2)-ENDMAT(didx,1)+1)
            tempdata = diodedata(1:endfreqidxs(didx),startidx:endidx)./...
                diodedata(1:endfreqidxs(didx),startidx-1);
            YDATA = [YDATA, reshape(tempdata,1,[])];
        end

        % Fit here
        sfunct=@(mua,mus) mdprepfun(...
            mua,mus,OPTS.nind,rhorange,DATAS.freqs(1:endfreqidxs(didx)),...
            startrhoidxs(didx),endidx);
        fitfunct=@(p,x)sfunct(p(1),p(2));
        
        fititer = fititer+1;
        [ops,~,~,OUTDATA.exits(fititer,didx)] = ...
            lsqcurvefit(fitfunct,[0.01,1.5],[],YDATA,[],[],options);
        OUTDATA.rmu(fititer,didx,:) = real(ops);
    end
    OUTDATA.exp{didx} = YDATA;
    optemp = median(squeeze(OUTDATA.rmu(:,didx,:)),1);
    OUTDATA.theory{didx} = sfunct(optemp(1),optemp(2));
end

% Fit scattering parameters
x = OPTS.laser_names(OPTS.usediodes);
pwrlaw = @(b,x) b(1).*(x.^-b(2));
for i = 1:length(OPTS.laser_names)
    y(i) = mean(OUTDATA.rmu(OUTDATA.exits(:,i)>0,i,2));
end
nrmrsd = @(b) norm(y(~isnan(y)) - pwrlaw(b,x(~isnan(y))));
OUTDATA.pwrfit = fminsearch(nrmrsd,[10000,1.4]);
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
lsqoptions = optimoptions('lsqcurvefit','Algorithm','Levenberg-Marquardt',...
    'FunctionTolerance',5e-12,'StepTolerance',1e-7,...
    'MaxFunctionEvaluations',60000,'MaxIterations',60000);
psoptions = optimoptions('patternsearch','FunctionTolerance',1e-14,...
    'MaxFunctionEvaluations',60000,'MaxIterations',60000,...
    'MeshTolerance',2*eps);

endfreqidxs = zeros(6,1);
startrhoidxs = endfreqidxs;
numtries = 5;
OUTDATA.exits = -1.*ones(numtries,length(OPTS.usediodes));
disp('lsqcurvefit Exit Flags:')

for didx = 1:length(OPTS.usediodes)
    trhorange = newrhos(didx,:);
    fititer = 0;
    
    % Round estimated endfreqs to existing frequencies
    endfreqidxs(didx) = find(DATAS.freqs>ENDMAT(didx,3),1,'first');
    startrhoidxs(didx) = find(OPTS.rhorange == ENDMAT(didx,1));
    
%     for endidx = (ENDMAT(didx,2)-4:ENDMAT(didx,2))-ENDMAT(didx,1)+1
    for jj = 1:numtries
        endidx = find(OPTS.rhorange==ENDMAT(didx,2))-numtries+jj;

        % For reshape operation to work as I would like, must do (elementwise)
        % transpose here.
        diodedata = DATAS.complex(:,1:endfreqidxs(didx),didx).';
        YDATA = [];
        for startidx = 1+startrhoidxs(didx):endidx
    %         for endidx = startidx+1:(ENDMAT(didx,2)-ENDMAT(didx,1)+1)
%             fprintf('%d%s%d\n',startidx,'...',endidx)
            tempdata = diodedata(:,startidx:endidx)./...
                diodedata(:,startidx-1);
            YDATA = [YDATA, reshape(tempdata,1,[])];
%             fprintf('%d\n',length(YDATA))
        end

        % Fit here
        sfunct=@(mua,mus) mdprepfun(...
            mua,mus,OPTS.nind,trhorange,DATAS.freqs(1:endfreqidxs(didx)),...
            startrhoidxs(didx),endidx);
        lsqfitfunct=@(p,x)sfunct(p(1),p(2));
        psfitfunct=@(p,x) sum(abs(sfunct(p(1),p(2))-YDATA));

        
        fititer = fititer+1;
        vals = [.02*rand(1), 2*rand(1)];
        tic
        [ops,fval(fititer,didx),OUTDATA.exits(fititer,didx),output] = ...
            patternsearch(psfitfunct,vals,[],[],[],[],[0,0],[.5,5],[],psoptions);
        toc
%         tic
%         [ops2,~,~,OUTDATA.exits2(fititer,didx)] = ...
%             lsqcurvefit(lsqfitfunct,vals,[],YDATA,[],[],lsqoptions);
%         toc
        OUTDATA.rmu(fititer,didx,:) = real(ops);
        format long, ops
        disp('PatternSearch Exit Flags:')
        OUTDATA.exits
        
%         figure;plot(YDATA,'.')
%         hold on;plot(sfunct(squeeze(OUTDATA.rmu(endidx-5,didx,1)),squeeze(OUTDATA.rmu(endidx-5,didx,2))),'x')
    end
    OUTDATA.exp{didx} = YDATA;
    if sum(OUTDATA.exits(:,didx) ~= 2) == 0
        OUTDATA.opfd(:,didx) = median(squeeze(OUTDATA.rmu(:,didx,:)));
    elseif sum(OUTDATA.exits(:,didx) == 2) == 1
        OUTDATA.opfd(:,didx) = OUTDATA.rmu(OUTDATA.exits(:,didx)==2,didx,:);
    else
        OUTDATA.opfd(:,didx) = median(squeeze(OUTDATA.rmu(OUTDATA.exits(:,didx)>0,didx,:)),1);
    end
    OUTDATA.theory{didx} = sfunct(OUTDATA.opfd(1,didx),OUTDATA.opfd(2,didx));
end
OUTDATA.endfidxs = endfreqidxs;

% Fit scattering parameters
x = OPTS.laser_names(OPTS.usediodes);
% pwrlaw = @(b,x) -b(1).*log(x)+b(2)+b(3).*x.^-4;
pwrlaw = @(b,x) b(1).*(x.^-b(2));
fdmusp = OUTDATA.opfd(2,:);

options = optimset('MaxFunEvals',10000,'Display','off');
nrmrsd = @(b) norm(fdmusp(~isnan(fdmusp)) - pwrlaw(b,x(~isnan(fdmusp))));
OUTDATA.pwrfit = fminsearch(nrmrsd,[4000,1.3],options);
% OUTDATA.pwrfit = fminsearch(nrmrsd,[1,10,1000],options);
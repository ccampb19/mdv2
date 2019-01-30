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
options = optimoptions('lsqcurvefit','Algorithm','Levenberg-Marquardt',...
    'FunctionTolerance',1e-12,'StepTolerance',1e-7,...
    'MaxFunctionEvaluations',60000,'MaxIterations',60000,...
    'Display','off');

endfreqidxs = zeros(6,1);
startrhoidxs = endfreqidxs;
numtries = 3;
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
        fitfunct=@(p,x)sfunct(p(1),p(2));
        
        fititer = fititer+1;
        [ops,~,~,OUTDATA.exits(fititer,didx)] = ...
            lsqcurvefit(fitfunct,[0.00522,1],[],YDATA,[],[],options);
        OUTDATA.rmu(fititer,didx,:) = real(ops);
        disp('lsqcurvefit Exit Flags:')
        OUTDATA.exits
        
%         figure;plot(YDATA,'.')
%         hold on;plot(sfunct(squeeze(OUTDATA.rmu(endidx-5,didx,1)),squeeze(OUTDATA.rmu(endidx-5,didx,2))),'x')
    end
    OUTDATA.exp{didx} = YDATA;
    if sum(OUTDATA.exits(:,didx) > 1) == 0
        opfd(:,didx) = median(squeeze(OUTDATA.rmu(:,didx,:)));
    elseif sum(OUTDATA.exits(:,didx)>1) == 1
        opfd(:,didx) = OUTDATA.rmu(OUTDATA.exits(:,didx)>0,didx,:);
    else
        opfd(:,didx) = mean(squeeze(OUTDATA.rmu(OUTDATA.exits(:,didx)>0,didx,:)),1);
    end
    OUTDATA.theory{didx} = sfunct(opfd(1,didx),opfd(2,didx));
end
OUTDATA.endfidxs = endfreqidxs;

% Fit scattering parameters
x = OPTS.laser_names(OPTS.usediodes);
pwrlaw = @(b,x) b(1).*(x.^-b(2));
fdmusp = opfd(2,:);

options = optimset('MaxFunEvals',1000,'Display','off');
nrmrsd = @(b) norm(fdmusp(~isnan(fdmusp)) - pwrlaw(b,x(~isnan(fdmusp))));
OUTDATA.pwrfit = fminsearch(nrmrsd,[8000,1.3],options);

% Do broadband fit?
% if OPTS.bb == 1
% %     Traditional broadband fit with scaling to FD MuA. I prefer to do
% %     multirho broadband with FD MuSP only, to confirm results are
% %     consistent.
%     muscat = pwrlaw(OUTDATA.pwrfit,DATAS.wv);
%     for didx = 1:length(OPTS.laser_names)
% %         fdwvidxs(i) = find(DATAS.wv>OPTS.laser_names(i),1,'first');
%         OUTDATA.fdmua(didx) = mean(OUTDATA.rmu(OUTDATA.exits(:,didx)>0,didx,1));
%     end
%     rref = interp1(DATAS.wv,DATAS.R,OPTS.laser_names);
%     rth = (Rtheory(OUTDATA.fdmua,fdmusp,OPTS.rhorange'-1.3,OPTS.nind))';
%     for ridx = 1:size(rth,2)
%         rscale(ridx) = rth(:,ridx)\rref(:,ridx);
%     end
%     % Trying a single scaling factor for now
%     rscaled = DATAS.R./median(rscale);
%     OUTDATA.bbmuas = zeros(size(rscaled));
%     disp('Calculating Broadband Reflectance...')
%     for ridx = 1:size(rscaled,2)
%         for widx = 1:size(rscaled,1)
% %                         OUTDATA.bbmuas(widx,ridx) = abs(fzero(@(mu) rscaled(widx,ridx) - ...
% 
%                         OUTDATA.bbmuas(widx,ridx) = abs(fzero(@(mu) DATAS.R(widx,ridx)./rscale(ridx) - ...
%                 abs(Rtheory(mu,muscat(widx),OPTS.rhorange(ridx)-1.3,OPTS.nind)),.01));
%         end
%     end
%     OUTDATA.wv = DATAS.wv;
%     disp('Done!')        
% end

if OPTS.bb == 1
    chopidxs = 425:1605;
    muscat = pwrlaw(OUTDATA.pwrfit,DATAS.wv);
    newrhos = OPTS.rhorange(1:end)-1;
    disp('Calculating Broadband Reflectance...')
    rchop = DATAS.R(chopidxs,1:end);
    muchop = muscat(chopidxs);
end
if OPTS.bb == 1
    % Fit for each wavelength
    for i = 1:size(rchop,1)
        for j = 1:5
        blah = rchop(i,2:end-j+1)./rchop(i,1:end-j);
        rfunct = @(mua,xdata) abs(Rtheory(mua,muchop(i),xdata(2:end),OPTS.nind))./...
                abs(Rtheory(mua,muchop(i),xdata(1:end-1),OPTS.nind));

        p1funct = @(mua,xdata) abs(p1seminfcompfit([mua,muchop(i)],0,0,OPTS.nind,xdata(2:end),0,0,1))./...
            abs(p1seminfcompfit([mua,muchop(i)],0,0,OPTS.nind,xdata(1:end-1),0,0,1));

        rfitteds(i,j) = abs(lsqcurvefit(rfunct,.01,newrhos(1:end-j+1),blah,[],[],options));
        p1fitteds(i,j) = abs(lsqcurvefit(p1funct,.01,newrhos(1:end-j+1),blah,[],[],options));
        end
    end
    OUTDATA.rfit = rfitteds;
    OUTDATA.p1fit = p1fitteds;
    OUTDATA.wv = DATAS.wv(chopidxs);
    disp('Done!')
end


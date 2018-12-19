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
    'FunctionTolerance',1e-9,'StepTolerance',1e-6,...
    'MaxFunctionEvaluations',50000,'MaxIterations',50000);

endfreqidxs = zeros(6,1);
startrhoidxs = endfreqidxs;
for didx = 1:length(OPTS.usediodes)
    trhorange = newrhos(didx,:);
    fititer = 0;
    for endidx = (ENDMAT(didx,2)-4:ENDMAT(didx,2))-ENDMAT(didx,1)+1

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
            mua,mus,OPTS.nind,trhorange,DATAS.freqs(1:endfreqidxs(didx)),...
            startrhoidxs(didx),endidx);
        fitfunct=@(p,x)sfunct(p(1),p(2));
        
        fititer = fititer+1;
        [ops,~,~,OUTDATA.exits(fititer,didx)] = ...
            lsqcurvefit(fitfunct,[0.022,1],[],YDATA,[],[],options);
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
    fdmusp(i) = mean(OUTDATA.rmu(OUTDATA.exits(:,i)>0,i,2));
end
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
    newrhos = OPTS.rhorange(1:end)-1.4;
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
                abs(Rtheory(mua,muchop(i),xdata(1:end-1),nchop(i)));

        p1funct = @(mua,xdata) abs(p1seminfcompfit([mua,muchop(i)],0,0,OPTS.nind,xdata(2:end),0,0,1))./...
            abs(p1seminfcompfit([mua,muchop(i)],0,0,OPTS.nind,xdata(1:end-1),0,0,1));

        rfitteds(i,j) = abs(lsqcurvefit(rfunct,.01,newrhos(1:end-j+1),blah,[],[],options));
        p1fitteds(i,j) = abs(lsqcurvefit(p1funct,.01,newrhos(1:end-j+1),blah,[],[],options));
        end
    end
    disp('Done!')
end
    
if OPTS.bb == 1
    % This experimental _if_ statement is an attempt at self-calibrated broadband.
    % Uniqueness of solutions has not been proven yet, so it may all be for
    % naught.
    reff = 2.1037*OPTS.nind^6-19.8048*OPTS.nind^5+76.8786*OPTS.nind^4-...
        156.9634*OPTS.nind^3+176.4549*OPTS.nind^2 -101.6004*OPTS.nind+22.9286;
    newrhos = OPTS.rhorange(1:end)'-1.3;
    chopidxs = 425:1180;
    options = optimset('MaxFunEvals',1e4,'MaxIter',1e4,'TolX',1e-9,'TolFun',1e-9,'Display','Iter');
    loptions = optimset('MaxFunEvals',1000,'Display','off');

    % Initial guess for preft & slope provided by FDPM for now
    x0 = [OUTDATA.pwrfit(1),OUTDATA.pwrfit(2)];
    for i = 1:2
        clear muass musparams
        for j = 1:4
        tic
        % Find solution for increasing numbers of points in a loop
        cchop = chopidxs(1:10^(3-i):end);
        wchop = DATAS.wv(cchop)';
        rchop = DATAS.R(cchop,1:end-j)';
        rrat = rchop(2:end,:)./rchop(1:end-1,:);
%         if i == 1
%             x0 = [x0,.01.*ones(size(cchop))];
%         end
        fitfun = @(x) mrhobb(x(1),x(2),rrat,newrhos(1:end-j),wchop,OPTS.nind,reff,loptions,0);
%         fitfun = @(x) sum(sum(sqrt((mrhobb(x(1),x(2),x(3:end),newrhos,...
%             wchop,OPTS.nind,reff)-rrat).^2)));
        [x,fval] = fminsearch(fitfun,x0,options);
        musparams(j,:) = x;
        [anss,muass(j,:)] = mrhobb(x(1),x(2),rrat,newrhos(1:end-j),wchop,OPTS.nind,reff,loptions,1);

        % Use solution as initial guess for next fitting iteration
        x0 = x;
        toc
        end
    end
    chopidxs2 = 425:1605;
    wchop2 = DATAS.wv(chopidxs2)';
    musps = x(1:2);
    muchop = pwrlaw(x(1:2),wchop2);
    rchop = DATAS.R(chopidxs2,1:14)';
    rrat = rchop(2:end,:)./rchop(1:end-1,:);
    for i = 1:6
        fdidxs(i) = find(DATAS.wv>OPTS.laser_names(i),1,'first');
    end
    %% is this the right way to do this?
    for i = 1:size(rchop,1)
        perfecttheory = Rtheory(mean(OUTDATA.rmu(:,:,1)),mean(OUTDATA.rmu(:,:,2)),...
            newrhos(i),OPTS.nind)';
        rscale(i) = perfecttheory\DATAS.R(fdidxs,i);
    end
    rrescale = DATAS.R(:,1:14).\rscale;
    for widx = 1:size(rrat,2)
%         blah = rrat(:,i);
        for ridx = 1:length(blah)
            sfunct = @(mua) sum(sqrt((abs(p1seminfcompfit([mua,muchop(i)],0,0,OPTS.nind,newrhos(2:end-j+1),0,0,1))./...
                abs(p1seminfcompfit([mua,muchop(i)],0,0,nchop(i),newrhos(1:end-j),0,0,1))...
                -blah).^2));
                    bbmuas(widx,ridx) = abs(fzero(@(mu) rrescale(widx,ridx) - ...
            abs(Rtheory(mu,muchop(widx),rho(ridx),n)),.01));
%             fcurve = @(mua,xdata) mrhobb(x0(1),x0(2),rrat,newrhos,wchop2,OPTS.nind,reff,loptions,1);
%             fitteds(i) = lsqcurvefit(fcurve,.005,newrhos,blah);
        end
    end
end



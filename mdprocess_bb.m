function OUTDATA = mdprocess_bb(OPTS,DATAS,OUTDATA)

pwrlaw = @(b,x) b(1).*(x.^-b(2));
chopidxs = 425:1605;
muscat = pwrlaw(OUTDATA.pwrfit,DATAS.wv);
newrhos = OPTS.rhorange(1:end)-1;
rchop = DATAS.R(chopidxs,1:end);
muchop = muscat(chopidxs);
wvchop = DATAS.wv(chopidxs);

% Traditional broadband fit with scaling to FD MuA. I prefer to do
% multirho broadband with FD MuS' only, to confirm results are
% consistent.
disp('Calculating Broadband Reflectance...')
for didx = 1:length(OPTS.laser_names)
%         fdwvidxs(i) = find(DATAS.wv>OPTS.laser_names(i),1,'first');
    OUTDATA.fdmua(didx) = mean(OUTDATA.rmu(OUTDATA.exits(:,didx)>0,didx,1));
end
rref = interp1(wvchop,rchop,OPTS.laser_names);
for i = 1:length(OUTDATA.fdmua)
    rth(i,:) = (p1seminfcompfit([OUTDATA.fdmua(i),OUTDATA.opfd(2,i)],0,0,OPTS.nind,newrhos',0,0,1));
end
for ridx = 1:size(rth,2)
    rscale(ridx) = rth(:,ridx)\rref(:,ridx);
end
% Scale each spectrum to theory
rscaled = rchop./rscale;
OUTDATA.rfit = zeros(size(rscaled));
for ridx = 1:size(rscaled,2)
    for widx = 1:size(rscaled,1)
%                         OUTDATA.bbmuas(widx,ridx) = abs(fzero(@(mu) rscaled(widx,ridx) - ...
%             OUTDATA.rfit(widx,ridx) = abs(fzero(@(mu) rchop(widx,ridx)./rscale(ridx) - ...
%                 abs(Rtheory(mu,muchop(widx),newrhos(ridx),OPTS.nind)),.01));

        OUTDATA.rfit(widx,ridx) = abs(fzero(@(mu) rscaled(widx,ridx) - ...
            abs(p1seminfcompfit([mu,muchop(widx)],0,0,OPTS.nind,newrhos(ridx),0,0,1)),.01));
    end
end
OUTDATA.wv = wvchop;
disp('Done!')


% Multirho broadband with FD MuS' only
disp('Calculating Broadband Reflectance...')
% options = optimset('MaxFunEvals',1000,'Display','off');
options = optimset('Display','off');
for i = 1:size(rchop,1)
    for j = 1:5
    blah = rchop(i,2:end-j+1)./rchop(i,1:end-j);
%         rfunct = @(mua,xdata) abs(Rtheory(mua,muchop(i),xdata(2:end),OPTS.nind))./...
%                 abs(Rtheory(mua,muchop(i),xdata(1:end-1),OPTS.nind));

    p1funct = @(mua,xdata) abs(p1seminfcompfit([mua,muchop(i)],0,0,OPTS.nind,xdata(2:end),0,0,1))./...
        abs(p1seminfcompfit([mua,muchop(i)],0,0,OPTS.nind,xdata(1:end-1),0,0,1));

%         rfitteds(i,j) = abs(lsqcurvefit(rfunct,.01,newrhos(1:end-j+1),blah,[],[],options));
    p1fitteds(i,j) = abs(lsqcurvefit(p1funct,.01,newrhos(1:end-j+1),blah,[],[],options));
    end
end
%     OUTDATA.rfit = rfitteds;
OUTDATA.p1fit = p1fitteds;
OUTDATA.wv = DATAS.wv(chopidxs);
disp('Done!')

% % This experimental _if_ statement is an attempt at self-calibrated broadband.
% % Uniqueness of solutions has not been proven yet, so it may all be for
% % naught.
% reff = 2.1037*OPTS.nind^6-19.8048*OPTS.nind^5+76.8786*OPTS.nind^4-...
%     156.9634*OPTS.nind^3+176.4549*OPTS.nind^2 -101.6004*OPTS.nind+22.9286;
% newrhos = OPTS.rhorange(1:end)'-1.3;
% chopidxs = 425:1180;
% options = optimset('MaxFunEvals',1e4,'MaxIter',1e4,'TolX',1e-9,'TolFun',1e-9,'Display','Iter');
% loptions = optimset('MaxFunEvals',1000,'Display','off');
% 
% % Initial guess for preft & slope provided by FDPM for now
% x0 = [OUTDATA.pwrfit(1),OUTDATA.pwrfit(2)];
% for i = 1:2
%     clear muass musparams
%     for j = 1:4
%         tic
%         % Find solution for increasing numbers of points in a loop
%         cchop = chopidxs(1:10^(3-i):end);
%         wchop = DATAS.wv(cchop)';
%         rchop = DATAS.R(cchop,1:end-j)';
%         rrat = rchop(2:end,:)./rchop(1:end-1,:);
% %         if i == 1
% %             x0 = [x0,.01.*ones(size(cchop))];
% %         end
%         fitfun = @(x) mrhobb(x(1),x(2),rrat,newrhos(1:end-j),wchop,OPTS.nind,reff,loptions,0);
% %         fitfun = @(x) sum(sum(sqrt((mrhobb(x(1),x(2),x(3:end),newrhos,...
% %             wchop,OPTS.nind,reff)-rrat).^2)));
%         [x,fval] = fminsearch(fitfun,x0,options);
%         musparams(j,:) = x;
%         [anss,muass(j,:)] = mrhobb(x(1),x(2),rrat,newrhos(1:end-j),wchop,OPTS.nind,reff,loptions,1);
% 
%         % Use solution as initial guess for next fitting iteration
%         x0 = x;
%         toc
%     end
% end
% chopidxs2 = 425:1605;
% wchop2 = DATAS.wv(chopidxs2)';
% musps = x(1:2);
% muchop = pwrlaw(x(1:2),wchop2);
% rchop = DATAS.R(chopidxs2,1:14)';
% rrat = rchop(2:end,:)./rchop(1:end-1,:);
% for i = 1:6
%     fdidxs(i) = find(DATAS.wv>OPTS.laser_names(i),1,'first');
% end
% %% is this the right way to do this?
% for i = 1:size(rchop,1)
%     perfecttheory = Rtheory(mean(OUTDATA.rmu(:,:,1)),mean(OUTDATA.rmu(:,:,2)),...
%         newrhos(i),OPTS.nind)';
%     rscale2(i) = perfecttheory\DATAS.R(fdidxs,i);
% end
% rrescale = DATAS.R(:,1:14).\rscale2;
% for widx = 1:size(rrat,2)
% %         blah = rrat(:,i);
%     for ridx = 1:length(blah)
%         sfunct = @(mua) sum(sqrt((abs(p1seminfcompfit([mua,muchop(i)],0,0,OPTS.nind,newrhos(2:end-j+1),0,0,1))./...
%             abs(p1seminfcompfit([mua,muchop(i)],0,0,nchop(i),newrhos(1:end-j),0,0,1))...
%             -blah).^2));
%                 bbmuas(widx,ridx) = abs(fzero(@(mu) rrescale(widx,ridx) - ...
%         abs(Rtheory(mu,muchop(widx),newrhos(ridx),n)),.01));
% %             fcurve = @(mua,xdata) mrhobb(x0(1),x0(2),rrat,newrhos,wchop2,OPTS.nind,reff,loptions,1);
% %             fitteds(i) = lsqcurvefit(fcurve,.005,newrhos,blah);
%     end
% end
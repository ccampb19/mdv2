function OUTDATA = mdprocess_bb(OPTS,DATAS,OUTDATA)

pwrlaw = @(b,x) b(1).*(x.^-b(2));
muscat = pwrlaw(OUTDATA.pwrfit,DATAS.wv);
% Adjust rho for fiber geometry (side-by-side 1mm fibers)
% rhocal = linspace(11,29,20);
newrhos = OPTS.bbrhorange-1.1;
rchop = DATAS.R;
muchop = muscat;
wvchop = DATAS.wv;

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
OUTDATA.rfit_sph = zeros(size(rscaled));
for ridx = 1:size(rscaled,2)
    for widx = 1:size(rscaled,1)
%                         OUTDATA.bbmuas(widx,ridx) = abs(fzero(@(mu) rscaled(widx,ridx) - ...
%             OUTDATA.rfit(widx,ridx) = abs(fzero(@(mu) rchop(widx,ridx)./rscale(ridx) - ...
%                 abs(Rtheory(mu,muchop(widx),newrhos(ridx),OPTS.nind)),.01));

        OUTDATA.rfit_sph(widx,ridx) = abs(fzero(@(mu) rscaled(widx,ridx) - ...
            abs(p1seminfcompfit([mu,muchop(widx)],0,0,OPTS.nind,newrhos(ridx),0,0,1)),.01));
    end
end
OUTDATA.rfit_sph(DATAS.bbdarkidxs) = 0;
OUTDATA.wv = wvchop;
disp('Done!')


% Multirho broadband with FD MuS' only
disp('Calculating Broadband Reflectance...')
options = optimset('MaxFunEvals',10000,'Display','off');
% psoptions = optimoptions('patternsearch','FunctionTolerance',1e-14,...
%     'MaxFunctionEvaluations',1800000,'MaxIterations',60000,...
%     'MeshTolerance',2*eps);
% options = optimset('Display','off');

% lsqoptions = optimoptions('lsqnonlin',...
%     'FunctionTolerance',5e-12,'StepTolerance',1e-7,...
%     'MaxFunctionEvaluations',60000,'MaxIterations',60000,...
%     'Display','off');

for i = 1:size(rchop,2)
    filteridxs{i} = DATAS.bbdarkcols(DATAS.bbdarkrows==i);
end
% rfitteds = zeros(size(rchop,1),size(rchop,2)-3);
rfitteds = zeros(size(rchop,1),DATAS.ncutoff-3);
p1fitteds = rfitteds;

    % Determine SDS pairs (uses unadjusted rhorange,outputs adjusted)
    % (For noise filtering, fixed later)
%     bswitch = 1;
%     srows = 1:DATAS.ncutoff;
%     nbrhos = newrhos(srows);

OUTDATA.bbpairs = mdrhopairs(newrhos(1:DATAS.ncutoff));
trhorange = newrhos(OUTDATA.bbpairs);
for i = 1:size(rchop,1)
    
    blindex = 1:DATAS.ncutoff;
    
    % Come back to reimplement noise filtering later.
    
%     srows = 1:DATAS.ncutoff;
%     bdata = rchop(i,srows);
%     for k = srows
%         if ~isempty(filteridxs{k})
%             if ismember(i,filteridxs{k})
%                 bdata(k) = nan;
%             end
%         end
%     end
%     if sum(~isnan(bdata))<4
%         rfitteds(i,j) = nan;
%         p1fitteds(i,j) = nan;
%     else
%         blindex = find(~isnan(bdata));
%         blah = bdata(blindex(2:end))./bdata(blindex(1:end-1));
%             debugRs{i,j} = blah;

%             rfunct = @(mua,xdata) abs(Rtheory(mua,muchop(i),xdata(2:end),OPTS.nind))./...
%                     abs(Rtheory(mua,muchop(i),xdata(1:end-1),OPTS.nind));

        zfun = @(mua,mus,xdata1,xdata2) ...
            Rtheory(mua,muchop(i),xdata1,OPTS.nind)./...
            Rtheory(mua,muchop(i),xdata2,OPTS.nind);                
        rvec = rchop(i,OUTDATA.bbpairs(:,1))./rchop(i,OUTDATA.bbpairs(:,2));
        zfunct = @(p,x) norm(zfun(p,muchop(i),trhorange(:,1),trhorange(:,2))-rvec);                

%             rfunct = @(p,xdata) p(2).*(Rtheory(p(1),muchop(i),xdata,OPTS.nind));
%             zfunct = @(pp) rfunct(pp,newrhos(blindex))-bdata(blindex);
%             temp = lsqnonlin(zfunct,[.01,5e8],[0,-inf],[10,inf],lsqoptions);
%             rfitteds(i,j) = temp(1);
%             OUTDATA.ayys(i,j) = temp(2);
        pfun = @(mua,xdata) ...
            p1seminfcompfit([mua,muchop(i)],0,0,OPTS.nind,xdata(:,1),0,0,1)./...
            p1seminfcompfit([mua,muchop(i)],0,0,OPTS.nind,xdata(:,2),0,0,1);
%         p1funct = @(mua,xdata) abs(p1seminfcompfit([mua,muchop(i)],0,0,OPTS.nind,xdata(2:end),0,0,1))./...
%             abs(p1seminfcompfit([mua,muchop(i)],0,0,OPTS.nind,xdata(1:end-1),0,0,1));

        [rfitteds(i),resnorms(i),~,exits(i)] = fminsearch(zfunct,.01);
%             rfitteds(i,j) = abs(rfitteds(i,j));

        p1fitteds(i) = abs(lsqcurvefit(pfun,.01,trhorange,rvec',[],[],options));

        %%%% DEBUG %%%%%%
        if ~mod(i,10)
%             if j==1 && (i==210 || i==877)
            figure(102)
            plot(rvec,'ko')
            hold on
            plot(zfun(rfitteds(i),muchop(i),trhorange(:,1),trhorange(:,2)));
            xlim([0, length(rvec)]);
            title(num2str(wvchop(i)))
            drawnow
            pause(.3)
            hold off

%                 figure(103)
%                 plot(abs(Rtheory(rfitteds(i,j),muchop(i),newrhos(blindex(2:end)),OPTS.nind)),'.');
%                 hold on
%                 plot(abs(Rtheory(rfitteds(i,j),muchop(i),newrhos(blindex(1:end-1)),OPTS.nind)),'.');
%                 title(num2str(wvchop(i)))

%         end
        %%%% /DEBUG %%%%%
    end
end

disp('Calculating Error')
OUTDATA.resnorms = resnorms;
OUTDATA.rerr = zeros(size(rfitteds,1),1);
OUTDATA.p1err = zeros(size(p1fitteds,1),1);
for i = 1:size(rfitteds,1)
    testd = rfitteds(i,:);
    OUTDATA.rerr(i) = std(testd(~isnan(testd)));
    testd = p1fitteds(i,:);
    OUTDATA.p1err(i) = std(testd(~isnan(testd)));
end
for i = 1:length(OPTS.laser_names)
    fbmismatchwvidxs(i) = min(find(wvchop<OPTS.laser_names(i),1,'last'),...
        find(wvchop>OPTS.laser_names(i),1,'first'));
end
OUTDATA.fbmismatch_r = mean(rfitteds(fbmismatchwvidxs,:),2)-OUTDATA.fdmua';
OUTDATA.fbmismatch_p1 = mean(p1fitteds(fbmismatchwvidxs,:),2)-OUTDATA.fdmua';

OUTDATA.rfit = rfitteds;
OUTDATA.p1fit = p1fitteds;
OUTDATA.wv = DATAS.wv;
disp('Done!')


% % This experimental _if_ statement is an attempt at self-calibrated broadband.
% % Uniqueness of solutions has not been proven yet, so it may all be for
% naught.
% % % % % reff = 2.1037*OPTS.nind^6-19.8048*OPTS.nind^5+76.8786*OPTS.nind^4-...
% % % % %     156.9634*OPTS.nind^3+176.4549*OPTS.nind^2 -101.6004*OPTS.nind+22.9286;
% % % % % wvidxs = 1:length(wvchop);
% % % % % options = optimset('MaxFunEvals',1e4,'MaxIter',1e4,'TolX',1e-9,'TolFun',1e-9,'Display','Iter');
% % % % % loptions = optimset('MaxFunEvals',1000,'Display','off');
% % % % % psoptions = optimoptions('patternsearch','FunctionTolerance',1e-14,...
% % % % %     'MaxFunctionEvaluations',1800000,'MaxIterations',60000,...
% % % % %     'MeshTolerance',2*eps);
% % % % % 
% % % % % 
% % % % % zfun = @(mua,mus,g) ...
% % % % %     Rtheoryg(mua,mus,newrhos(1:DATAS.ncutoff),OPTS.nind,g,0,reff);
% % % % % rvec = mdprepbb(rchop(1,1:DATAS.ncutoff));
% % % % % zfunct = @(p,x) norm(mdprepbb(zfun(p(1),p(2),p(3)))-rvec);
% % % % % OUTDATA.gfit(1,:) = patternsearch(zfunct,[.01,10*muchop(5),0.8],...
% % % % %         [],[],[],[],[1e-6,10*(muchop(5)-.2),.4],[.07,10*(muchop(5)+.2),.93],[],psoptions);
% % % % %     
% % % % % for i = 1:size(rchop,1)/5-1
% % % % %     wwidx = 5*i+5;
% % % % %     if ~mod(i,50)
% % % % %         disp(i);
% % % % %     end
% % % % %     rvec = mdprepbb(rchop(wwidx,1:DATAS.ncutoff));
% % % % %     zfunct = @(p,x) norm(mdprepbb(zfun(p(1),p(2),p(3)))-rvec);
% % % % %     OUTDATA.gfit(i+1,:) = patternsearch(zfunct,[OUTDATA.gfit(i,1),OUTDATA.gfit(i,2),OUTDATA.gfit(i,3)],...
% % % % %         [],[],[],[],[1e-6,10*(muchop(wwidx)-.2),.4],[.07,10*(muchop(wwidx)+.2),.93],[],psoptions);
% % % % % end

% Initial guess for preft & slope provided by FDPM
% for i = 1:2
%     clear muass musparams
%         tic
%         % Find solution for increasing numbers of points in a loop
%     cchop = wvidxs(1:10^(3-i):end);
%     wchop = DATAS.wv(cchop)';
%     rchop = DATAS.R(cchop,1:end)';
%     rrat = rchop(2:end,:)./rchop(1:end-1,:);
% %         if i == 1
% %             x0 = [x0,.01.*ones(size(cchop))];
% %         end
%     fitfun = @(x) mrhobb(x(1),x(2),rrat,newrhos(1:end-j),wchop,OPTS.nind,reff,loptions,0);
% %         fitfun = @(x) sum(sum(sqrt((mrhobb(x(1),x(2),x(3:end),newrhos,...
% %             wchop,OPTS.nind,reff)-rrat).^2)));
%     [x,fval] = fminsearch(fitfun,x0,options);
%     musparams(j,:) = x;
%     [anss,muass(j,:)] = mrhobb(x(1),x(2),rrat,newrhos(1:end-j),wchop,OPTS.nind,reff,loptions,1);
% 
%     % Use solution as initial guess for next fitting iteration
%     x0 = x;
%     toc
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
%% is this the right way to do this?
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
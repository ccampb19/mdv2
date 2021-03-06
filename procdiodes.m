function OUTDATA = procdiodes(ENDMAT,DATAS,OPTS,newrhos,psoptions)

OUTDATA.exits = -1.*ones(1,length(OPTS.usediodes));
endfreqidxs = zeros(6,1);
startrhoidxs = endfreqidxs;

for didx = 1:length(OPTS.usediodes)

    fititer = 0;

    % Round estimated endfreqs to existing frequencies
    endfreqidxs(didx) = find(DATAS.freqs>=ENDMAT(didx,3),1,'first');
    startrhoidxs(didx) = find(OPTS.rhorange == ENDMAT(didx,1));

    %     for endidx = (ENDMAT(didx,2)-4:ENDMAT(didx,2))-ENDMAT(didx,1)+1
    endidx = find(OPTS.rhorange==ENDMAT(didx,2));
    
    % Determine SDS pairs (uses unadjusted rhorange,outputs adjusted)
    OUTDATA.fdpairs = mdrhopairs(OPTS.rhorange(startrhoidxs(didx):endidx));
    trhorange = reshape(newrhos(didx,OUTDATA.fdpairs),[],2);

    % For reshape operation to work as I would like, must do (elementwise)
    % transpose here.
    diodedata = DATAS.complex(:,1:endfreqidxs(didx),didx).';
    YDATA = [];
    for yidx = 1:size(OUTDATA.fdpairs,1)
        tempdata = diodedata(:,OUTDATA.fdpairs(yidx,1))./...
            diodedata(:,OUTDATA.fdpairs(yidx,2));
        YDATA = [YDATA, reshape(tempdata,1,[])];
    end

    % Fit here
    sfunct=@(mua,mus) mdprepfun(...
        mua,mus,OPTS.nind,trhorange,DATAS.freqs(1:endfreqidxs(didx)));
%     lsqfitfunct=@(p,x) sfunct(abs(p(1)),abs(p(2)))-YDATA;
    psfitfunct=@(p,x) sum(abs(sfunct(p(1),p(2))-YDATA)./abs(YDATA));


    fititer = fititer+1;
    vals = [.02*rand(1), 2*rand(1)];
    tic
    [ops,fval(fititer,didx),OUTDATA.exits(fititer,didx),output] = ...
        patternsearch(psfitfunct,vals,[],[],[],[],[0,0],[.5,5],[],psoptions);
    toc
    
%     tic
%     [ops2,~,~,OUTDATA.exits2(fititer,didx)] = ...
%         lsqcurvefit(lsqfitfunct,vals,[],YDATA,[],[],lsqoptions);
%     toc
    
    OUTDATA.rmu(fititer,didx,:) = real(ops);
    format long, ops
    disp('PatternSearch Exit Flags:')
    OUTDATA.exits

    OUTDATA.exp{didx} = YDATA;
    OUTDATA.opfd(:,didx) = ops;
%     if sum(OUTDATA.exits(:,didx) ~= 2) == 0
%         OUTDATA.opfd(:,didx) = median(squeeze(OUTDATA.rmu(:,didx,:)));
%     elseif sum(OUTDATA.exits(:,didx) == 2) == 1
%         OUTDATA.opfd(:,didx) = OUTDATA.rmu(OUTDATA.exits(:,didx)==2,didx,:);
%     else
%         OUTDATA.opfd(:,didx) = median(squeeze(OUTDATA.rmu(OUTDATA.exits(:,didx)>0,didx,:)),1);
%     end
    OUTDATA.theory{didx} = sfunct(OUTDATA.opfd(1,didx),OUTDATA.opfd(2,didx));
end

OUTDATA.endfidxs = endfreqidxs;

end
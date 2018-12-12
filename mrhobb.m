function [msum,muas] = mrhobb(A,B,RRAT,RHO,LAMBDA,NIND,REFF,LOPT,FLAG)
    persistent mua0
    if isempty(mua0)
        mua0 = .01.*ones(size(LAMBDA));
    end

    muas = 0;
    musp = A.*LAMBDA.^(-B);
    fitteds = zeros(size(musp));
    resn = fitteds;
    for i = 1:length(musp)
        rratwv = RRAT(:,i);
        fcurve = @(mua,xdata) Rtheory(mua,musp(i),xdata(2:end),NIND(i),0,REFF(i))./...
            Rtheory(mua,musp(i),xdata(1:end-1),NIND(i),0,REFF(i));
        if ~FLAG
            [~,resn(i)] = lsqcurvefit(fcurve,.005,RHO,rratwv,[],[],LOPT);
        else
            [muas(i),resn(i)] = lsqcurvefit(fcurve,.005,RHO,rratwv,[],[],LOPT);
        end
    end
    msum = sum(resn);        
end  
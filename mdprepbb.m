function outvec = mdprepbb(RDATA,RHOS)

outvec = [];
for tdidx = 1:size(RHOS,1)
    outvec = [outvec, RDATA(RHOS(tdidx,1))/RDATA(RHOS(tdidx,2))];
end
    


function MODELDATA = mdprepfun(MUA,MUSP,NIND,RHOS,FREQS)

p = [MUA,MUSP];
MODELDATA=[];
for tdidx = 1:size(RHOS,1)
    tempdata = p1seminfcompfit(p,FREQS,0,NIND,RHOS(tdidx,1),0,0,1)./...
        p1seminfcompfit(p,FREQS,0,NIND,RHOS(tdidx,2),0,0,1);
    MODELDATA = [MODELDATA, reshape(tempdata,1,[])];
end
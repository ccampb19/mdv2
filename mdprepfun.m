function MODELDATA = mdprepfun(MUA,MUSP,NIND,RHOS,FREQS,STARTRHO,ENDRHO)

p = [MUA,MUSP];

nrho = length(RHOS);

MODELDATA=[];

for startidx = (STARTRHO+1):(ENDRHO-1)
    tempdata = p1seminfcompfit(p,FREQS,0,NIND,RHOS(startidx:ENDRHO),0,0,1)./...
        p1seminfcompfit(p,FREQS,0,NIND,RHOS(startidx-1),0,0,1);
    MODELDATA = [MODELDATA, reshape(tempdata,1,[])];
end
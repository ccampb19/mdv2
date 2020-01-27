function rp = mdrhopairs(rhorange)

mxsep = rhorange(end)-rhorange(1);
rp = [];
if mxsep < 3
    for i = 1:length(rhorange)-1
        for j = i+1:length(rhorange)
            rp = [rp;j,i];
        end
    end
elseif mxsep < 5
    for i = 1:length(rhorange)-1
        trs = find(rhorange >= rhorange(i)+3);
        for j = 1:length(trs)
            rp = [rp;trs(j),i];
        end
    end
else
    for i = 1:length(rhorange)-1
        trs = find(rhorange >= rhorange(i)+5);
        for j = 1:length(trs)
            rp = [rp;trs(j),i];
        end
    end    
end
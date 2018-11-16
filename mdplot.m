function mdplot(OPTS,OUTPUT)
%MDPLOT makes plots, extremely simplified for now
%   OPTS contains all necessary data to load a mutlidistance dataset.
%   ENDMAT, 3 columns, has start/end rhos and end frequencies
%   EXITS contains exitflags indicating degree of fitting success

diodenum = length(OPTS.usediodes);

figure
for i = 1:diodenum
    subplot(diodenum,2,2*i-1)
    hold on
    plot(real(OUTPUT.exp{i}),'.')
    plot(real(OUTPUT.theory{i}))
    ylabel('Real')
    title(num2str(OPTS.laser_names(OPTS.usediodes(i))))
    
    subplot(diodenum,2,2*i)
    hold on
    plot(imag(OUTPUT.exp{i}),'.')
    plot(imag(OUTPUT.theory{i}))
    ylabel('Imag.')
    title(num2str(OPTS.laser_names(i)))
end

figure
subplot(2,1,1)
hold on
for i = 1:5
plot(OPTS.laser_names(OUTPUT.exits(i,:)>0),OUTPUT.rmu(i,(OUTPUT.exits(i,:)>0),1),'x')
end

ylabel('\mu_A (mm^-^1)')
title('Calculated OPs')

subplot(2,1,2)
hold on
for i = 1:5
    plot(OPTS.laser_names(OUTPUT.exits(i,:)>0),OUTPUT.rmu(i,(OUTPUT.exits(i,:)>0),2),'x')
end
plot(600:1000,OUTPUT.pwrfit(1).*(600:1000).^-OUTPUT.pwrfit(2),'k--')
xlabel('Wavelength (nm)')
ylabel('\mu_S'' (mm^-^1)')
    
% for i = 1:6
%     muas = OUTPUT.rmu(:,i,1);
%     muss = OUTPUT.rmu(:,i,2);
%     try
%         cov(i,1) = 100.*std(muas(OUTPUT.exits(:,i)>0))./mean(muas(OUTPUT.exits(:,i)>0));
%         cov(i,2) = 100.*std(muss(OUTPUT.exits(:,i)>0))./mean(muss(OUTPUT.exits(:,i)>0));
%     catch
%         cov(i,:) = 0;
%     end
%     cov2(i,1) = 100.*std(muas)./mean(muas);
%     cov2(i,2) = 100.*std(muss)./mean(muss);
% end
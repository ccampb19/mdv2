function phantom = mdplot(OPTS,OUTPUT)
%MDPLOT makes plots, extremely simplified for now
%   OPTS contains all necessary data to load a mutlidistance dataset.
%   ENDMAT, 3 columns, has start/end rhos and end frequencies
%   EXITS contains exitflags indicating degree of fitting success

diodenum = length(OPTS.usediodes);
colors = varycolor_rainbow(length(OPTS.rhorange));

figure
for i = 1:diodenum
    tempe = reshape(OUTPUT.exp{i},OUTPUT.endfidxs(i),[]);
    tempt = reshape(OUTPUT.theory{i},OUTPUT.endfidxs(i),[]);
    
    subplot(diodenum,2,2*i-1)
    hold on
    ha = plot(abs(tempe(1:10,:))','.');
    hb = plot(abs(tempt(1:10,:))');
    set(hb,{'Color'},num2cell(colors(1:10,:),2));
    ylabel('Amp.')
    title(num2str(OPTS.laser_names(OPTS.usediodes(i))))
    
    subplot(diodenum,2,2*i)
    hold on
    ha = plot(angle(tempe(1:10,:)'),'.','Color',colors(i,:));
    set(ha,{'Color'},num2cell(colors(1:10,:),2));
    hb = plot(angle(tempt(1:10,:)'),'Color',colors(i,:));
    set(hb,{'Color'},num2cell(colors(1:10,:),2));
    ylabel('Phase')
    title(num2str(OPTS.laser_names(i)))
    
end
suptitle('Phase/Amp Ratios & Fits')

figure(1000)
subplot(2,1,1)
hold on
for i = 1:size(OUTPUT.rmu,1)
    plot(OPTS.laser_names(OUTPUT.exits(i,:)>0),OUTPUT.rmu(i,(OUTPUT.exits(i,:)>0),1),'x')
    plot(OPTS.laser_names(OUTPUT.exits(i,:)==0),OUTPUT.rmu(i,(OUTPUT.exits(i,:)==0),1),'o')
end
xlim([600,1000])

ylabel('\mu_A (mm^-^1)')
title('Calculated OPs')

subplot(2,1,2)
hold on
for i = 1:size(OUTPUT.rmu,1)
    plot(OPTS.laser_names(OUTPUT.exits(i,:)>0),OUTPUT.rmu(i,(OUTPUT.exits(i,:)>0),2),'x')
    plot(OPTS.laser_names(OUTPUT.exits(i,:)==0),OUTPUT.rmu(i,(OUTPUT.exits(i,:)==0),2),'o')
end
plot(600:1000,OUTPUT.pwrfit(1).*(600:1000).^-OUTPUT.pwrfit(2),'k--')
xlabel('Wavelength (nm)')
ylabel('\mu_S'' (mm^-^1)')

if OPTS.bb == 1
%     figure
    figure(1000)
    subplot(2,1,1)
    hold on
    for i = 1:size(OUTPUT.rmu,1)
        plot(OPTS.laser_names(OUTPUT.exits(i,:)>0),OUTPUT.rmu(i,(OUTPUT.exits(i,:)>0),1),'kx')
        plot(OPTS.laser_names(OUTPUT.exits(i,:)==0),OUTPUT.rmu(i,(OUTPUT.exits(i,:)==0),1),'ro')
    end
    
    plot(OUTPUT.wv,OUTPUT.rfit,':');
    plot(OUTPUT.wv,OUTPUT.p1fit);
    
    xlabel('\lambda (nm)')
    ylabel('\mu_A (mm^-^1)')
end

figure
for i = 1:length(OPTS.usediodes)
    tempexpmat = reshape(OUTPUT.exp{i},OUTPUT.endfidxs(i),[]);
    tempthymat = reshape(OUTPUT.theory{i},size(tempexpmat));
    subplot(2,floor(length(OPTS.usediodes)/2),i)
    title([num2str(OPTS.laser_names(OPTS.usediodes((i)))) ' nm'])
    hold on
    for j = 1:length(OPTS.rhorange)
%     for j = 1:size(tempexpmat,2)
        plot(tempexpmat(:,j),'.','Color',colors(j,:))
        plot(tempthymat(:,j),'.','Color',colors(j,:),'MarkerSize',1)
    end        
    xlabel('Re')
    ylabel('Im')
end
suptitle('FD Complex Fits')

for i = 1:length(OPTS.usediodes)
    muas(i) = mean(OUTPUT.rmu(OUTPUT.exits(:,i)>0,i,1));
    muss(i) = mean(OUTPUT.rmu(OUTPUT.exits(:,i)>0,i,2));
end
muas(isnan(muas)) = median(OUTPUT.rmu(:,isnan(muas),1));
muss(isnan(muss)) = median(OUTPUT.rmu(:,isnan(muss),2));
phantom = [OPTS.laser_names;muas;muss]';

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
function [phantom, bbp, bap] = mdplot(OPTS,OUTPUT)
%MDPLOT makes plots, extremely simplified for now
%   OPTS contains all necessary data to load a mutlidistance dataset.
%   ENDMAT, 3 columns, has start/end rhos and end frequencies
%   EXITS contains exitflags indicating degree of fitting success

diodenum = length(OPTS.usediodes);
numlines = min(10,length(OPTS.rhorange));
colors = varycolor_rainbow(numlines);

% figure
% for i = 1:diodenum
%     tempe = reshape(OUTPUT.exp{i},OUTPUT.endfidxs(i),[]);
%     tempt = reshape(OUTPUT.theory{i},OUTPUT.endfidxs(i),[]);
%     
%     subplot(diodenum,2,2*i-1)
%     hold on
%     ha = plot(abs(tempe(1:numlines,:))','.');
%     hb = plot(abs(tempt(1:numlines,:))');
%     set(hb,{'Color'},num2cell(colors(1:numlines,:),2));
%     ylabel('Amp.')
%     title(num2str(OPTS.laser_names(OPTS.usediodes(i))))
%     
%     subplot(diodenum,2,2*i)
%     hold on
%     ha = plot(angle(tempe(1:numlines,:)'),'.','Color',colors(i,:));
%     set(ha,{'Color'},num2cell(colors(1:numlines,:),2));
%     hb = plot(angle(tempt(1:numlines,:)'),'Color',colors(i,:));
%     set(hb,{'Color'},num2cell(colors(1:numlines,:),2));
%     ylabel('Phase')
%     title(num2str(OPTS.laser_names(i)))
%     
% end
% sgtitle('Phase/Amp Ratios & Fits')

fignum = 1000+randi(1000);
figure(fignum)
subplot(2,1,1)
hold on
for i = 1:size(OUTPUT.rmu,1)
    plot(OPTS.laser_names(OUTPUT.exits(i,:)==2),OUTPUT.rmu(i,(OUTPUT.exits(i,:)==2),1),'kx')
    plot(OPTS.laser_names(OUTPUT.exits(i,:)~=2),OUTPUT.rmu(i,(OUTPUT.exits(i,:)~=2),1),'ro')
end
xlim([600,1000])

ylabel('\mu_A (mm^-^1)')
title([OPTS.filenameprototype ' Calculated OPs'],'Interpreter','none')

subplot(2,1,2)
hold on
for i = 1:size(OUTPUT.rmu,1)
    plot(OPTS.laser_names(OUTPUT.exits(i,:)==2),OUTPUT.rmu(i,(OUTPUT.exits(i,:)==2),2),'kx')
    plot(OPTS.laser_names(OUTPUT.exits(i,:)~=2),OUTPUT.rmu(i,(OUTPUT.exits(i,:)~=2),2),'ro')
end
% plot(600:1000,-OUTPUT.pwrfit(1).*log(600:1000)+OUTPUT.pwrfit(2),'k--')
plot(600:1000,OUTPUT.pwrfit(1).*(600:1000).^-OUTPUT.pwrfit(2),'k--')
xlabel('Wavelength (nm)')
ylabel('\mu_S'' (mm^-^1)')

if OPTS.bb == 1
%     figure
    figure(fignum)
    subplot(2,1,1)
    hold on
    
    ha = plot(OUTPUT.wv,OUTPUT.rfit,'--');
    set(ha,{'Color'},num2cell(varycolor_rainbow(size(OUTPUT.rfit,2)),2));
%     plot(OUTPUT.wv,OUTPUT.rfit_sph);
    hb = plot(OUTPUT.wv,OUTPUT.p1fit);
    set(hb,{'Color'},num2cell(varycolor_rainbow(size(OUTPUT.p1fit,2)),2));
    
    for i = 1:size(OUTPUT.rmu,1)
        plot(OPTS.laser_names(OUTPUT.exits(i,:)==2),OUTPUT.rmu(i,(OUTPUT.exits(i,:)==2),1),'kx')
        plot(OPTS.laser_names(OUTPUT.exits(i,:)~=2),OUTPUT.rmu(i,(OUTPUT.exits(i,:)~=2),1),'ro')
    end
    
    xlabel('\lambda (nm)')
    ylabel('\mu_A (mm^-^1)')
end

figure
for i = 1:length(OPTS.usediodes)
    tempexpmat = reshape(OUTPUT.exp{i},OUTPUT.endfidxs(i),[]);
    tempthymat = reshape(OUTPUT.theory{i},size(tempexpmat));
    subplot(2,ceil(length(OPTS.usediodes)/2),i)
    title([num2str(OPTS.laser_names(OPTS.usediodes((i)))) ' nm'])
    hold on
    for j = 1:min([length(OPTS.rhorange),size(tempexpmat,2),10])
%     for j = 1:size(tempexpmat,2)
        if ~isempty(tempexpmat)
            plot(tempexpmat(:,j),'.','Color',colors(j,:))
            plot(tempthymat(:,j),'.','Color',colors(j,:),'MarkerSize',1)
        end
    end        
    xlabel('Re')
    ylabel('Im')
end
sgtitle('FD Complex Fits')

for i = 1:length(OPTS.usediodes)
    muas(i) = mean(OUTPUT.rmu(OUTPUT.exits(:,i)>0,i,1));
    muss(i) = mean(OUTPUT.rmu(OUTPUT.exits(:,i)>0,i,2));
end
muas(isnan(muas)) = median(OUTPUT.rmu(:,isnan(muas),1));
muss(isnan(muss)) = median(OUTPUT.rmu(:,isnan(muss),2));
phantom = [OPTS.laser_names;muas;muss]';

if OPTS.bb == 1
    for i = 1:length(OUTPUT.p1fit)
        tvec = OUTPUT.p1fit(i,:);
    % for i = 1:length(OUTPUT.p1fit)
    %     tvec = OUTPUT.p1fit(i,:);
        if sum(tvec) == 0
            means(i) = 0;
            ameans(i) = 0;
        elseif length(tvec) < 3
            means(i) = mean(tvec(tvec~=0));
            ameans(i) = means(i);
        else
            means(i) = mean(tvec(tvec~=0));
            ameans(i) = mean(tvec(end-2:end));
        end
    end
    bbp = [OUTPUT.wv(:),means(:),OUTPUT.pwrfit(1).*OUTPUT.wv(:).^-OUTPUT.pwrfit(2)];
    bap = [OUTPUT.wv(:),ameans(:),OUTPUT.pwrfit(1).*OUTPUT.wv(:).^-OUTPUT.pwrfit(2)];
else
    bbp = 0;
    bap = 0;
end
end

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
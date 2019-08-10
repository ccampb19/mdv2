function [fit,concs,wvrange] = mdchromfit(phantom,chromfile)

temp = importdata(chromfile);
chromwv = temp.data(:,1);
chromdata = temp.data(:,2:end);
chromnames=temp.colheaders(2:end)';

wvrange = phantom(:,1);
A = interp1(chromwv,chromdata,wvrange);
pidxs = ~isnan(phantom(:,2)) & ~isnan(A(:,1));
muawv = phantom(pidxs,1);
%Don't fit past 910nm
testwv = find(muawv>910,1,'first');
if isempty(testwv)
    testwv = length(mua);
end
pidxs(testwv:end) = 0;

mua = phantom(pidxs,2);
muawv = phantom(pidxs,1);

concs = lsqnonneg(A(pidxs,:),mua);
fit = chromdata*concs;

figure
plot(muawv,mua,chromwv,fit)
xlabel('\lambda (nm)')
ylabel('\mu_A (mm^-^1)')
legend('Recov. \mu_A','Fit')

disp(['StO2: '  num2str(100*concs(1)./(concs(1)+concs(2))) '%']);
disp(['Total Hb: ' num2str(concs(1)+concs(2))]);
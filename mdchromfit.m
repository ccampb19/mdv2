function [fit,concs,wvrange] = mdchromfit(phantom,chromfile)

temp = importdata(chromfile);
chromwv = temp.data(:,1);
chromdata = temp.data(:,2:end);
chromnames=temp.colheaders(2:end)';

wvrange = phantom(~isnan(phantom(:,2)),1);
mua = phantom(~isnan(phantom(:,2)),2);
A = interp1(chromwv,chromdata,wvrange);

concs = lsqnonneg(A,mua);
fit = chromdata*concs;

figure
plot(wvrange,mua,chromwv,fit)
xlabel('\lambda (nm)')
ylabel('\mu_A (mm^-^1)')
legend('Recov. \mu_A','Fit')
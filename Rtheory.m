function R = Rtheory(mua, musp, rho, n, boundary_opt, Reff)

% 3/06 AEC updated with proper Reff and new integral calc done by Sophie

%return the steday-state reflectance - semi-infinite case
% mua=absorption coeff
% musp= reduce scattering coeff
% rho=source-detector separation
% n=refractive index

% Reff=-1.440/n^2+0.710/n+0.668+0.0636*n;			%assumes air-tissue interface
% Reff = 2.1037.*n.^6-19.8048.*n.^5+76.8786.*n.^4-156.9634.*n.^3+176.4549.*n.^2 -101.6004.*n+22.9286;
 if nargin<5 
     boundary_opt=0; 
 end
 if nargin<6
     if n==1.4
         Reff=0.493;
     elseif n==1.33
         Reff=0.431;
     else
         if boundary_opt==1
             Reff=Ref_n_lookup_v2(n);   %use integrals
         else    %polynomial fit 6 order by Sophie and AEC
             Reff = 2.1037.*n.^6-19.8048.*n.^5+76.8786.*n.^4-156.9634.*n.^3+176.4549.*n.^2 -101.6004.*n+22.9286;
         end
     end
 end


mutp=mua+musp;				% reduced att. coefficient

alfa = 1-0.8.*(musp+mua)/(musp.*1.6+mua);
D=1./(3.*(musp+alfa.*mua));			% diffusion constant, low albedo
% D=1./(3.*(musp+mua));       % diffusion constant
% D=1./(3*musp);      % Least accurate, probably?
mueff=(mua./D).^0.5;		% effective att.
% zo=1./musp;
zo=1./(mutp);				% distance surface-isotropic source
zb=2*D*(1+Reff)/(1-Reff); %distance  surface-extrapolated boundary
r1=(zo.^2+rho.^2).^0.5;  %
r2=((2*zb+zo).^2+rho.^2).^0.5; %

fluence=((4*pi*D).^-1).*(exp(-mueff.*r1)./r1 -exp(-mueff.*r2)./r2); %fluence
flux=((4*pi).^-1).*(zo.*(mueff+r1.^-1).*(exp(-mueff.*r1)./r1.^2) +((zo+2*zb).*(mueff+r2.^-1).*exp(-mueff.*r2)./r2.^2));

fr1=0.7857 * (n.^3) - 4.3259 * (n.^2) + 8.26405 * n - 4.71306;
fr2=-0.02043 * (n.^3)- 0.38418 * (n.^2) + 2.01132 * n - 1.62198;
R=(1-fr1)/4*fluence-(fr2-1)/2*flux;  % see Kienle and Patterson, JOSA A 14(1), 246-254,1997
% R = fluence;

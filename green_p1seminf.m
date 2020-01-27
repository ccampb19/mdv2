%========================================================================================
%SEMI-INFINITE PHASE AND AMPLITUDE CALCULATOR FOR GREEN P1 MODEL, v1.0 
%========================================================================================
% Returns complex value of P1 seminfinite reflectance.  See Ola for more details.
%
% USAGE
%    [Y, VER] = GREEN_P1SEMINF(P,F,FX,NIND,RHO1,RHO2,WT,REIM_FLAG,OPT1);  
%
% Note on Units:    Accepts MHz for frequency, 1/mm for optical properties, and mm
%   for distances
%
%INPUT:	accepts an array of Nx1 frequencies (one r) 
%               P is the input optical properties in the format  P=[mua,mus], both in 1/mm 
%				F is the modulation frequency in MHz (meant to be the primary independent variable 
%				FX is an option to depermine what should be fit.  FX=0 fits both mua and mus. 
%                   If FX==-1, fits mua only, and if FX =+1, fits mus only. <DEACTIVATED> . 
%				NIND is the sample index of refraction
%				RHO1 is the source-detector separation on surface in mm
%               RHO2 is the remainder of the distances.  Set RHO2=0 for single distance fit (default).  
%                   Only use nonzero RHO2 if multi-distance fit is required.  Make RHO2 a vector if there are
%                   more than 2 distances <NOT OPERATIONAL>
%               WT is a weight to fit the frequencies; set to 0 if none are
%                   desired.
%               REIM_FLAG is binary marker: set to 0 to fit Real and Imiginary components, 
%                   set to 1 to fit Phase and Amplitude  (default) 
%               BOUND_OPT indicates the choice of boundary condition based upon  Haskell, R. C., L. O. Svaasand, T. Tsong-Tseh, 
%                   F. Ti-Chen, M. S. McAdams and B. J. Tromberg (1994). "Boundary conditions for the diffusion equation in radiative transfer." 
%                   Journal of the Optical Society of America A (Optics, Image Science and Vision) 11(10): 2727-41.
%                   Set to 0 to use precalculated values, 1 to calculate directly (slower ...)
%
%OUTPUT:	
%           Y returns a 2Nx1 matrix for, where 
%				fa		row 1..N    is amplitude P1 approximation PDW
%				fb		row N+1..2N is phase P1 approximation PDW (radians)
%           VER returns the version number


function [y, VER] = green_p1seminf(p,f,fx,n,rho,~,wt,reim_flag,~,~,p1flag)
VER = '1.0';

assert(size(f,2) == 1, 'Frequency vector must be of size nx1');
assert(size(rho,1) == 1, 'Rho vector must be of size(1xn');

if nargin<11
    p1flag = 1;
end

c=3e8*1000; %mm/s
cn=c/n;
ft = 1e6.*f;

mua = real(p(1));
musp = real(p(2));

mutr=mua+musp;
D=1/(3*mutr); 
zp=1/mutr;

A=-0.13755 * (n.^3) + 4.3390 * (n.^2) - 4.90466 * n + 1.68960;
zb = 2 / (3 * mutr) * A;

fr1=0.7857 * (n.^3) - 4.3259 * (n.^2) + 8.26405 * n - 4.71306;
fr2=-0.02043 * (n.^3)- 0.38418 * (n.^2) + 2.01132 * n - 1.62198;

y = zeros(length(ft),length(rho));
for fidx = 1:length(ft)
    
    if(p1flag == 0)
        % %     SDA k-vector
        k = sqrt(((mua * cn + 1i * ft(fidx) * 2 * pi) /(cn * D)));
    elseif(p1flag == 1)    
        % % %     p1 k-vector
        fbc = ft(fidx)*2*pi/cn;
        alpha = 3*fbc*D;
        k = sqrt((mua-fbc.*alpha-1i.*(fbc+mua*alpha))./D);

    else
        error('no p1flag')
    end

   zr1= -zp;
   zr2= zp+2*zb;

   r11=sqrt(rho.^2+ zp.^2);    %calculating the r of the source from the origin (0,0)
   r12=sqrt(rho.^2 +(zp+2*zb).^2);  %calculating the r of the probed point from the origin (0,0)

   [fluence,flux]=TemporalFrequencyGreenFunction5(D,k,r11,r12,zr1,zr2,0);

    y(fidx,:) = ((1 - fr1)/4).*fluence + ((fr2 - 1)/2).*flux;
    

end

if reim_flag == 0
    y = [abs(y),angle(y)];
end

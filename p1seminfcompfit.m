%========================================================================================
%SEMI-INFINITE PHASE AND AMPLITUDE CALCULATOR FOR P1 MODEL, v 2.5
%========================================================================================
% Returns phase (radians) and amplitude of P1 seminfinite Photon density
% waves (PDW).  See Cerussi, A. and B. Tromberg (2003), Photon Migration Spectroscopy: 
% Frequency Domain Techniques. Biomedical Photonics Handbook. T. Vo-Dinh,
% CRC Press: 22.1- 22.17 for more details.
%
% USAGE
%    [Y, VER] = P1SEMINF(P,F,FX,NIND,RHO1,RHO2,WT,REIM_FLAG,OPT1);  
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
%
%NOTES
%
%10/07 AEC poly fit used to make Reff calculation a lookup.  Error is less
%               than 0.002 in Reff calculation
%
%5/07 AEC Modified with extra argument to use full calc of reflectance from
%               Haskell as calculated by Sophie.  Emperical boundary
%               condition scrapped.
%
%1/00 AEC	Uses emperical value for reflection coefficient from 
%				R.A.J. Groenhuis et al. Appl. Opt. 22 2463 (1983).		
%
%1/00 AEC	Modified the Kr and Ki functions to eliminate the branch point at high
%				frequencies by entering Josh's Thesis P1 value for k (page 109) and letting
%				MATLAB take the real and imaginary parts.				


function y = p1seminfcompfit(p,f,fx,nind,rho1,rho2, wt, reim_flag,boundary_opt, reff)

%*************************************************************************
%PRELIMINARY MATTERS
%*************************************************************************

%set version number
%ver=2.5;

%argument check
if nargin <9
    boundary_opt=0; % Use polynomial approximation for haskell bounday condition
end

if nargin<8  %equivalent to phase/amp fit
   reim_flag = 0;
end

if nargin<7  %equivalent to no weights
   wt = ones(length(f)*2,1);
end

if nargin<6  %equivalent to one distance
   rho2=0;
end

%Other checking of inputs
if fx==-1	%fix mua then fit mus
    mua = real(p(1));			   %fixed, and taken in mm-1
    mus = real(p(2));				%fitted, taken in mm-1
elseif fx==-1					%fix mus then fit mua
    mua = real(p(1));			   %fixed, and taken in mm-1
    mus = real(p(2));				%fitted, taken in mm-1
else
    mua = real(p(1));			   %fixed, and taken in mm-1
    mus = real(p(2));				%fitted, taken in mm-1
end


%*************************************************************************
%BASIC CALCULATIONS
%*************************************************************************

%definitions 
c = 2.99792458e11/nind;				% now in mm/s
mutr = mua+mus;
ltr = 1/mutr;
D=1/3*ltr;			%diffusion coefficient, uses the mua for kicks
I=sqrt(-1);
fbc = 1.0e6*2*pi*f./c;			% now in MHZ for omega
alpha=3*fbc*D;		%Josh definition, such that alpha is 2pi*Tcoll/Tmod (page 109).

%boundary conditions
if nargin<10
    if nind==1.4
        reff=0.493;
    elseif nind==1.33
        reff=0.431;
    else
        if boundary_opt==1
            reff=Ref_n_lookup_v2(nind);   %use integrals
        else    %polynomial fit 6 order by Sophie and AEC
            reff = 2.1037.*nind.^6-19.8048.*nind.^5+76.8786.*nind.^4-156.9634.*nind.^3+176.4549.*nind.^2 -101.6004.*nind+22.9286;
        end
    end
end

%calcuate true distances   
zb = 2/3*(1+reff)/(1-reff)*ltr;									%extrapolated boundary

r01 = sqrt(ltr*ltr+rho1.*rho1);  	            %s-d separation for source, dist 1   
rb1 = sqrt((2*zb+ltr)*(2*zb+ltr)+rho1.*rho1);	%s-d separation for image, dist 1

%===kvector=====================================================================================
k_josh=sqrt((mua-fbc.*alpha-I.*(fbc+mua*alpha))./D);	%complete P1 wavevector
kr=abs(real(k_josh));	
ki=abs(imag(k_josh));		

%==photon density wave===========================================================================
er01 = exp(-kr.*r01);					%exponentials in form e(-kr)/r
erb1 = exp(-kr.*rb1);
Re1 = (+(er01./r01).*cos(ki.*r01) - (erb1./rb1).*cos(ki.*rb1))./D;
Im1 = (+(er01./r01).*sin(ki.*r01) - (erb1./rb1).*sin(ki.*rb1))./D;  


%*************************************************************************
%OUTPUT
%*************************************************************************
if(rho2==0)		%flag for single distance = real and imaginary
	if(reim_flag == 0)
        error('Complex fit requires fdpm.opt.imagfit = 1')
% 		fa=sqrt(Re1.^2+Im1.^2);
% 		fb=unwrap(atan2(Im1,Re1));		%in radians
	else
		fa=Re1;
		fb=Im1;
    end
else					%in phase and amplitude
   r02 = sqrt(ltr*ltr+rho2.*rho2);  	            %s-d separation for source, dist 2   
   rb2 = sqrt((2*zb+ltr)*(2*zb+ltr)+rho2.*rho2);	%s-d separation for image, dist 2
   er02 = exp(-kr.*r02);		
   erb2 = exp(-kr.*rb2);
   Re2 = (+(er02./r02).*cos(ki.*r02) - (erb2./rb2).*cos(ki.*rb2))./D;
   Im2 = (+(er02./r02).*sin(ki.*r02) - (erb2./rb2).*sin(ki.*rb2))./D;  
   
   if(reim_flag==0)
       error('Complex fit requires fdpm.opt.imagfit = 1')
% 	   fa=sqrt(Re1.^2+Im1.^2)./sqrt(Re2.^2+Im2.^2);
% 	   fb=unwrap(atan2(Im1,Re1))-unwrap(atan2(Im2,Re2));		%in radians
   else
	   fa=Re2;
	   fb=Im2;
   end
end

%final outputs
% y=[fa,fb];		%final return
y=fa+1i.*fb;

 if(wt~=0)
     y = y.*wt;  %Weight for the curvefitting
 end
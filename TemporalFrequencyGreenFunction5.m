function [F,J]=TemporalFrequencyGreenFunction5(D,k,r1,r2,zr1,zr2,del)
%Complex Temporal Frequency Point Source Image Greens Function fluence
Fluence_s1=(exp(-(k) .* r1+1i*del) ./ (4.0 * pi * D*r1));
Fluence_i1=(exp(-(k) .* r2+1i*del) ./ (4.0 * pi * D*r2));
F=Fluence_s1-Fluence_i1;
    
%Complex Temporal Frequency Point Source Image Greens Function flux
Flux_s1=(zr1./((r1.^2).*4.*pi)) .* (k + (1 ./ r1)) .* exp(-k.* r1+1i*del);
Flux_i1=(zr2./((r2.^2).*4.*pi)) .* (k + (1 ./ r2)) .* exp(-k.* r2+1i*del);
J=Flux_s1-Flux_i1;


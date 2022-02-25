%hm=phi/N0.*exp(j*2*pi/N0*n); 
%hpm=F0*phip/N0.*exp(j*2*pi/N0*n);;
%hppm=F0^2*phipp/N0.*exp(j*2*pi/N0*n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,p,ap,pp,app,ppp]=...
statestimator(hm,hpm,hppm,s)
%hm,hpm,hppm impulse response of the 
%FIR filter and first two derivatives
%sh samples of the filtered signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%a,p       amplitude and phase
%ap,pp     first derivatives of a and p
%app,ppp   second derivvatives of a and p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Amplitude and Phase estimation
k=(0:Nc*N0-1)'; p=conv(s, hm,'same');
ah=abs(p); ph=angle(p);   % rotating phase
phc=angle(p.*exp(-2j*pi/N0*k));% antirotated 
a=2*ah; p=phc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Frequency estimation
pp=conv(s,hpm,'same');
ppc=pp.*exp(-j*ph);        
ahp=real(ppc);     
php=imag(ppc)./ah;
ap=ahp; pp=php/(2*pi); %frequency deviation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ROCOF Estimation
ppp=conv(s,hppm,'same');
pppc=ppp.*exp(-j*ph); 
ahpp=real(pppc)+ah.*php.^2;  
phpp=( imag(pppc)- 2*ahp.*php )./ah;
app=ahpp; ppp=phpp/(2*pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
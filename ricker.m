function p=ricker(dt,fdom,tlength)
%function [wavelet,tw]=ricker(dt,fdom,tlength)
% creates a Ricker wavelet

if(nargin<3)
   tlength=127.*dt;
 end
 if(nargin<2)
   fdom=15.; 
 end
% create a time vector
  nt=round(tlength/dt)+1;
  %tmin=-dt*round(nt/2);
 % tmin=-180*dt;   
 tmin=-200*dt;
 tw=tmin+dt*(0:nt-1)';
% create the wavelet
  pf=pi^2*fdom^2;
  wavelet=(1-2.*pf*tw.^2).*exp(-pf*tw.^2);

% normalize
% generate a refenence sinusoid at the dominant frequency
%wavelet=wavenorm(wavelet,tw,2);
p=wavelet;
plot(p);

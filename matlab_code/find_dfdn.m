function dfdn = find_dfdn(f,nn,L)
% ASEN 5327 Class Notes Fall 2010 
N=nn-1;
n=[0:(N/2)-1 -(N/2):-1];
%L=1; %show fuction is periodic on wavelength of L
kn=2*pi*n./L;
kn((N/2)+1)=0;

%find Fourier components
f_tilde=fft(f,N);
dfdn=ifft(i*kn.*f_tilde,N);
dfdn(nn)=dfdn(1);

end


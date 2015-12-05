function  d2fdn2 = find_d2fdn2(f,nn,L) 
%ASEN5327 Class notes Fall 2010 
N=nn-1;
n=[0:(N/2)-1 -(N/2):-1];
kn=2*pi*n./L;
kn((N/2)+1)=0; 
f_tilde=fft(f,N);
d2fdn2=ifft(-(kn.^2).*f_tilde,N);
d2fdn2(nn)=d2fdn2(1);
end



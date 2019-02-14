a = abs(FFT2D((exp(1000*1i*phi))));
b = 20*log(a.*2);
g = improfile(b,[1,1023],[512,512]);
plot(x,g)
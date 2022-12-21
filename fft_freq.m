function freq=fft_freq(n,dx)
if mod(n,2)==0
freq=[0:(n/2-1),-n/2:-1]/(n*dx);
else
freq=[0:(n-1)/2,-(n-1)/2:-1]/(n*dx); 
end

end
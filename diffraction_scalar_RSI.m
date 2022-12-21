function  Uxyz=diffraction_scalar_RSI(wavelength,n,Uz0,xygird,zlist,OneDorTwoD)

% Diffraction simulation according to the first kind Rayleigh Sommerfield Integral
% 公式来自Handbook of optical system 的 sclar diffraction章节
% The simulation is based on Fast Fourier Methods
% ---------------------------------------------------------------------------------
% wavelength: input wave wavelength
% n: refractive index
% Uz0: The input wave's complex amplitude
% xy_grid: The x,y,z grid functions: {xgrid,ygrid,zgrid}
% z_list : The z points we want to simulate
% OneDorTwoD: 1D or 2D cases

zlist=[zlist];
zlist_1d=zlist(:);

if OneDorTwoD=="2D"
    xgrid=xygird{1};
    ygrid=xygird{2};
    dx=xgrid(1,2)-xgrid(1,1);
    dy=ygrid(2,1)-ygrid(1,1);
    [numy,numx]=size(xgrid);
    numz=size(zlist_1d,1);

    Uxyz=zeros(numy,numx,numz);

    Uz0_FFT=fft2(Uz0);

    [num_size_y,num_size_x]=size(Uz0_FFT);
    freq_x=fft_freq(num_size_x,dx);
    freq_y=fft_freq(num_size_y,dy);
    [freq_grid_x,freq_grid_y]=meshgrid(freq_x,freq_y);

    k=1/wavelength*n;
    kz=sqrt(k^2-freq_grid_x.^2-freq_grid_y.^2);
   
    for l=1:numz
        z=zlist_1d(l);
        PI=exp(2*pi*1i*kz*z);
        Uxyz(:,:,l)=ifft2(Uz0_FFT.*PI);
    end

end

if OneDorTwoD=="1D"
    xgrid=xygird{1};
    dx=xgrid(2)-xgrid(1);
    numx=size(xgrid(:),1);

    numz=size(zlist_1d,1);

    Uxyz=zeros(numx,numz);

    Uz0_FFT=fft(Uz0);

    freq_x=fft_freq(numx,dx);

    k=1/wavelength*n;
    kz=sqrt(k.^2-freq_x.^2-0);
   
    for l=1:numz
        z=zlist_1d(l);
        PI=exp(2*pi*1i*kz*z);
        Uxyz(:,l)=ifft(Uz0_FFT(:).*PI(:));
    end

end
end
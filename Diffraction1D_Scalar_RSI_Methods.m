%% 1D case
% -------------------------------------------------------------------------
nm=1e-9;
width_x=500*nm;
width_y=500*nm;
wavelength=100*nm;
n=1;

num_f=301;
k0=1/wavelength*n;
deltaf=k0*2/num_f;

dx=0.01*wavelength;
Tmax=1/deltaf;


x_supp=-Tmax/2:dx:Tmax/2;
numx=size(x_supp(:),1);


numz=601;
z_supp=linspace(0,30*wavelength,numz);

[xz_zgrid,xz_xgrid]=meshgrid(z_supp/wavelength,x_supp/wavelength);


% Rectangle
Uz0=ones(numx,1);
Uz0(abs(x_supp)>width_x/2)=0;


% The FFT Methods
xy_grid={x_supp};
tic
Uxz_FFT=diffraction_scalar_RSI(wavelength,n,Uz0,xy_grid,z_supp,"1D");
t2=toc;

%%
figure()
pcolor(xz_zgrid,xz_xgrid,abs(Uxz_FFT).^2);
colormap jet;
shading interp;
ylim([-10,10]);
colorbar();
title('Diffraction 1D results')
xlabel('x/\lambda');
ylabel('y/\lambda');
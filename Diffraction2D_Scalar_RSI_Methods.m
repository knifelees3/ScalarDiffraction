nm=1e-9;
n=1;

width_x=500*nm;
width_y=500*nm;
wavelength=100*nm;
num_f=100;
k=1/wavelength*n;
fmax=2*k/num_f;
deltaf=fmax/num_f;

% Tmax=0.001/deltaf;
Tmax=20*wavelength;
dx=0.01*wavelength;
dy=dx;
x_supp=-Tmax/2:dx:Tmax/2;
y_supp=-Tmax/2:dy:Tmax/2;

numx=size(x_supp(:),1);
numy=numx;
numz=1;
z_supp=[1000]*nm;

[xgrid,ygrid,zgrid]=meshgrid(x_supp,y_supp,z_supp);
xyz_grid={xgrid,ygrid,zgrid};

% Rectangle
Uz0=ones(numy,numx,numz);
Uz0(abs(xgrid)<width_x/2 &abs(ygrid)<width_y/2)=0;


% The FFT Methods
xy_grid={xgrid,ygrid};
tic
Uxyz_FFT=diffraction_scalar_RSI(wavelength,n,Uz0(:,:,1),xy_grid,z_supp,"2D");
t2=toc;

[xxgrid,yygrid]=meshgrid(x_supp,y_supp);

figure()
sgtitle(['Z/\lambda=10'])
subplot(121)
pcolor(xxgrid,yygrid,Uz0);
colormap jet; shading interp;
colorbar;
xlabel('x/\lambda');
ylabel('y/\lambda');

subplot(122)
pcolor(xxgrid,yygrid,abs(Uxyz_FFT));
colormap jet; shading interp;
colorbar;
xlabel('x/\lambda');
ylabel('y/\lambda');
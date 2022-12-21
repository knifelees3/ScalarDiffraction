nm=1e-9;
numx=301;
x_supp=linspace(-1000,1000,numx)*nm;
dx=x_supp(2)-x_supp(1);

numy=301;
y_supp=linspace(-1000,1000,numy)*nm;
dy=y_supp(2)-y_supp(1);

width_x=500*nm;
width_y=500*nm;
wavelength=100*nm;
n=1;
k=2*pi/wavelength*n;

numz=101;
z_supp=linspace(0,wavelength*10,numz);

[xgrid,ygrid,zgrid]=meshgrid(x_supp,y_supp,z_supp);
xyz_grid={xgrid,ygrid,zgrid};

% Rectangle
Uz0=ones(numy,numx,numz);
Uz0(abs(xgrid)>width_x/2)=0;
num_care=101;
x_care=linspace(-width_x,width_x,num_care);
y_care=[0];

xy_supp={x_care,y_care};

Uz0(abs(xgrid)>width_x/2)=0;
theta=pi/3;

Uz0_1=Uz0;
Uz0_2=Uz0.*exp(-1i*k*sin(theta).*xgrid);
Uz0_3=Uz0.*exp(1i*k*sin(theta).*xgrid);

xyz_grid={xgrid,ygrid,zgrid};
xy_supp={x_care,y_care};
Uxyz_1=diffraction_scalar(wavelength,n,Uz0_1,xyz_grid,xy_supp,"RSI");
Uxyz_2=diffraction_scalar(wavelength,n,Uz0_2,xyz_grid,xy_supp,"RSI");
Uxyz_3=diffraction_scalar(wavelength,n,Uz0_3,xyz_grid,xy_supp,"RSI");

[xz_xgrid,xz_zgrid]=meshgrid(x_care,z_supp);
Uxyz_1_2D=squeeze(abs(Uxyz_1(1,:,:)));
Uxyz_2_2D=squeeze(abs(Uxyz_2(1,:,:)));
Uxyz_3_2D=squeeze(abs(Uxyz_3(1,:,:)));
Uxyz_1_2D(:,1)=ones(num_care,1);
Uxyz_2_2D(:,1)=ones(num_care,1);
Uxyz_3_2D(:,1)=ones(num_care,1);
Uxyz_1_2D(abs(xz_xgrid)>width_x)=0;
Uxyz_2_2D(abs(xz_xgrid)>width_x)=0;
Uxyz_3_2D(abs(xz_xgrid)>width_x)=0;
%%
figure()
sgtitle(['D/\lambda=',num2str(width_x/wavelength,2)]);
subplot(311)
pcolor(xz_xgrid/wavelength,xz_zgrid/wavelength,Uxyz_1_2D);
colormap hot; shading interp;
colorbar;
% caxis([0,0.02]);
xlabel('x/\lambda');
ylabel('y/\lambda');
title(['\theta=',num2str(0,2)]);
subplot(312)
pcolor(xz_xgrid/wavelength,xz_zgrid/wavelength,Uxyz_2_2D);
colormap hot; shading interp;
colorbar;
% caxis([0,0.02]);
xlabel('x/\lambda');
ylabel('y/\lambda');
title(['\theta=-\pi/3']);

subplot(313)
pcolor(xz_xgrid/wavelength,xz_zgrid/wavelength,Uxyz_3_2D);
colormap hot; shading interp;
colorbar;
% caxis([0,0.02]);
xlabel('x/\lambda');
ylabel('y/\lambda');
title(['\theta=\pi/3']);
set(gcf,'position',[100 100 1200 800])
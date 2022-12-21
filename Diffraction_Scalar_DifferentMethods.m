nm=1e-9;
numx=201;
x_supp=linspace(-1000,1000,numx)*nm;
dx=x_supp(2)-x_supp(1);

numy=301;
y_supp=linspace(-1000,1000,numy)*nm;
dy=y_supp(2)-y_supp(1);

width_x=500*nm;
width_y=500*nm;
wavelength=100*nm;

numz=1;
z_supp=[1000]*nm;

[xgrid,ygrid,zgrid]=meshgrid(x_supp,y_supp,z_supp);
xyz_grid={xgrid,ygrid,zgrid};

% Rectangle
Uz0=ones(numy,numx,numz);
Uz0(abs(xgrid)>width_x/2)=0;
Uz0(abs(ygrid)>width_y/2)=0;

x_care=linspace(-width_x,width_x,21);
y_care=linspace(-width_y,width_y,21);

xy_supp={x_care,y_care};
Uxyz_RSI=diffraction_scalar(wavelength,n,Uz0,xyz_grid,xy_supp,"RSI");
Uxyz_Fresnel=diffraction_scalar(wavelength,n,Uz0,xyz_grid,xy_supp,"Fresnel");
Uxyz_Fraunhofer=diffraction_scalar(wavelength,n,Uz0,xyz_grid,xy_supp,"Fraunhofer");

[xy_xgrid,xy_ygrid]=meshgrid(x_care,y_care);
%%
lz=1;
figure()
sgtitle(['D/\lambda=',num2str(width_y/wavelength,2),'  z/\lambda=',num2str(z_supp(lz)/wavelength)]);
subplot(131);
pcolor(xy_xgrid/wavelength,xy_ygrid/wavelength,abs(Uxyz_RSI(:,:,1)));
colormap hot; shading interp;
colorbar;
caxis([0,1]);
xlabel('x/\lambda');
ylabel('y/\lambda');
title('RSI')
subplot(132);
pcolor(xy_xgrid/wavelength,xy_ygrid/wavelength,abs(Uxyz_Fresnel(:,:,1)));
colormap hot; shading interp;
colorbar;
caxis([0,1]);
xlabel('x/\lambda');
ylabel('y/\lambda');
title('Fresnel')
subplot(133);
pcolor(xy_xgrid/wavelength,xy_ygrid/wavelength,abs(Uxyz_Fraunhofer(:,:,1)));
colormap hot; shading interp;
colorbar;
caxis([0,1]);
xlabel('x/\lambda');
ylabel('y/\lambda');
title('Fraunhofer')
set(gcf,'units','normalized','position',[0.1 0.1 0.6 0.3]);

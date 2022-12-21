function Uxyz=diffraction_scalar(wavelength,n,Uz0,xyz_grid,xy_supp,method)

k0=2*pi/wavelength;
k=k0*n;

% 2D Case ： 公式来自Handbook of optical system 的 sclar diffraction章节
% ---------------------------------------------------------------------------------
% if OneDorTwoD=="2D"

% wavelength: input wave wavelength
% n: refractive index
% Uz0: The input wave's complex amplitude
% xyz_grid: The x,y,z grid functions: {xgrid,ygrid,zgrid}
% xy_supp : the x and y points that we care {xgrid_want,ygrid_want}
% methods: RSI: First Rayleigh–Sommerfeld diffraction
%          Fresnel: Fresnel diffraction
%          Fraunhofer: Fraunhofer Diffraction


xgrid=xyz_grid{1};
ygrid=xyz_grid{2};
zgrid=xyz_grid{3};
[~,~,numz]=size(xgrid);
x_supp=xy_supp{1};
y_supp=xy_supp{2};
numx=size(x_supp(:),1);
numy=size(y_supp(:),1);

x_supp_int=transpose(squeeze(xgrid(1,:,1)));
y_supp_int=squeeze(ygrid(:,1,1));

Uzz_RSI=zeros(numy,numx,numz);
Uzz_Fres=zeros(numy,numx,numz);
Uzz_Fraun=zeros(numy,numx,numz);
count=0;
if method=="RSI" 
    f = waitbar(0,"Simulating the diffraction with RSI methods");
    for m=1:numy
        for l=1:numx
            count=count+1;
            waitbar(count/numx/numy,f);
            Dis_RSI=sqrt(abs(xgrid-x_supp(l)).^2+abs(ygrid-y_supp(m)).^2+abs(zgrid).^2);
            costheta=zgrid./Dis_RSI;
            tobeint_RSI=Uz0.*(1i*k-1./Dis_RSI).*exp(1i*k*Dis_RSI)./Dis_RSI.*costheta/2/pi;
            Uzz_RSI(m,l,:)=trapz(y_supp_int,trapz(x_supp_int,tobeint_RSI,2));
        end
    end
    Uxyz=Uzz_RSI;
end    
    
if method=="Fresnel"
    f = waitbar(0,"Simulating the diffraction with Fresnel methods");
    for m=1:numy
        for l=1:numx
            count=count+1;
            waitbar(count/numx/numy,f);
            Dis_Fres=((xgrid-x_supp(l)).^2+(ygrid-y_supp(m)).^2)./zgrid/2;
            tobeint_Fres=Uz0.*exp(1i*k*Dis_Fres)*1i*k.*exp(1i*k*zgrid)./zgrid/2/pi;
            Uzz_Fres(m,l,:)=trapz(y_supp_int,trapz(x_supp_int,tobeint_Fres,2));
        end
    end
    Uxyz=Uzz_Fres;  
end

if method=="Fraunhofer"
    f = waitbar(0,"Simulating the diffraction with Fraunhofer methods");
    for m=1:numy
        for l=1:numx
                count=count+1;
                waitbar(count/numx/numy,f);
                Dis_Fraun=-(xgrid*x_supp(l)+ygrid*y_supp(m))./zgrid;
                tobeint_Fraun=Uz0.*exp(1i*k*Dis_Fraun).*1i.*exp(1i*k*zgrid)/wavelength./abs(zgrid).*exp(1i*k*(x_supp(l).^2+y_supp(m).^2)/2./zgrid);
                Uzz_Fraun(m,l,:)=trapz(y_supp_int,trapz(x_supp_int,tobeint_Fraun,2));
        end
    end
    Uxyz=Uzz_Fraun;
end

% if zgrid(1,1,1)==0
% Uxyz(:,:,1)=Uz0(:,:,1);
% end

end
% 一维计算下面这种方法是错误的，要计算一维需要认为将另外一个反向变成无穷大
% if OneDorTwoD=="1D"
%     xgrid=xyz_grid{1};
%     zgrid=xyz_grid{2};
%     dx=xgrid(1,2)-xgrid(1,1);
%     [numz,~]=size(xgrid);
%     x_supp=xy_supp{1};
%     z_supp=zgrid(:,1);
%     numx=size(x_supp(:),1);
%     x_supp_int=squeeze(xgrid(1,:,1));
%     Uzz_RSI=zeros(numz,numx);
%     Uzz_Fres=zeros(numz,numx);
%     Uzz_Fraun=zeros(numz,numx);
%     
% if method=="RSI"
%     f = waitbar(0,"Simulating the diffraction with RSI methods 1D case");
%     for m=1:numx
%             waitbar(m/numx,f);
%             Dis_RSI=sqrt(abs(xgrid-x_supp(m)).^2+abs(zgrid).^2);
%             costheta=zgrid./Dis_RSI;
%             tobeint_RSI=(ones(numz,1)*Uz0).*(1i*k-1./Dis_RSI).*exp(1i*k*Dis_RSI)./Dis_RSI.*costheta/2/pi;
%             Uzz_RSI(:,m)=trapz(x_supp_int,tobeint_RSI,2)*dx;
%     end
%     Uxyz=Uzz_RSI;
% end
% 
% 
% if method=="Fresnel"
%     f = waitbar(0,"Simulating the diffraction with Fresnel methods 1D case");
%     for m=1:numx
%             waitbar(m/numx,f);
%             Dis_Fres=(xgrid-x_supp(m)).^2./zgrid/2;
%             tobeint_Fres=(ones(numz,1)*Uz0).*exp(1i*k*Dis_Fres);
%             Uzz_Fres(:,m)=trapz(x_supp_int,tobeint_Fres,2)*1i*k.*exp(1i*k*z_supp(:))./z_supp(:)/2/pi*dx;
%     end
%    Uxyz=Uzz_Fres;
% end
% 
%     
% if method=="Fraunhofer"
%     f = waitbar(0,"Simulating the diffraction with Fraunhofer methods 1D case");
%     for m=1:numx
%             waitbar(m/numx,f);
%             Dis_Fraun=-xgrid*x_supp(m)./zgrid;
%             tobeint_Fraun=(ones(numz,1)*Uz0).*exp(1i*k*Dis_Fraun);
%             Uzz_Fraun(:,m)=trapz(x_supp_int,tobeint_Fraun,2).*1i.*exp(1i*k*z_supp(:))/...
%                 wavelength./abs(z_supp(:)).*exp(1i*k*x_supp(m).^2/2./z_supp(:))*dx;
%     end
%     Uxyz=Uzz_Fraun;
% end
%    
% end




% end

function [] = Phiplot(filename,prec)

% Phiplot(filename,prec)
% Potential contour plots of filename='... .dat'. prec is an integer
% indicating the meshing coarsness. prec=1: large file. prec=2 or prec=4:
% smaller file.

figure
readoutput;

[R,T]=meshgrid(rcc,tcc);

XR=R.*T;
YR=R.*sqrt(1-T.*T);

Dtheta1=acos(0.25*(3+tcc(2)));
Dtheta2=acos(tcc(2));
Dr=(rcc(nrused)-1)./(nrused-1);

[r,the]=meshgrid(1:Dr*prec:rcc(nrused)+Dr/nrused,0:pi/round(pi/(prec*Dtheta1)):pi);
xr=r.*cos(the);
yr=r.*sin(the);

% modify potential values at theta=0 and pi in order to account for the fact
% that the first cell center is not at thete=0 or pi. Use 2d order forward
% derivative =0 on axis

A=[1 1 1 ; 0 Dtheta1 Dtheta2 ; 0 Dtheta1^2 Dtheta2^2];
B=[0;1;0];
C=A^-1*B;
for k=1:nrused
    phi(k,1)=-(C(2)*phi(k,1)+C(3)*phi(k,2))/C(1);
    phi(k,nthused)=-(C(2)*phi(k,nthused)+C(3)*phi(k,nthused-1))/C(1);
end
phir=griddata(XR,YR,phi',xr,yr);
pcolor(xr,yr,phir);
shading interp

hold all


%[r,the]=meshgrid(2:Dr*16*prec:rcc(nrused)+Dr/nrused, ...
%    pi/round(pi/(prec*Dtheta1)):pi/round(pi/(prec*Dtheta2)):pi-pi/round(pi/(prec*Dtheta1)));
%xr=r.*cos(the);
%yr=r.*sin(the);
%vx=griddata(XR,YR,(vrsum.*T-vtsum.*sqrt(1-T.^2))./(psum+0.001),xr,yr);
%vy=griddata(XR,YR,(vtsum.*T+vrsum.*sqrt(1-T.^2))./(psum+0.001),xr,yr);

%quiver(xr,yr,vx,vy,1,'k')


[xx,yy]=meshgrid(1:3,0:1/63.:1);

[XX,YY]=meshgrid(1:3,power((0:1/127.:1),2));
comap=griddata(xx,yy,jet,XX,YY);
caxis([-0.5 0]);
colormap(comap);

axis equal
axis([-rcc(nrused) rcc(nrused) 0 rcc(nrused)]);

colorbar('NorthOutside')
hold on
%X=-1:0.01:1;
%Y=sqrt(1-power(X,2));
%plot(rmax*X,rmax*Y,'k','LineWidth',2)

[C,h]=contour(xr,yr,phir,[-0.15 -0.085 -0.05],'k','LineWidth',1);
text_handle=clabel(C,h);
set(text_handle,'FontSize',16);



%[C,h]=contour(xr,yr,phir,[-0.3],'k','LineWidth',1);
%text_handle=clabel(C,h);
%set(text_handle,'FontSize',16);

ylabel('\rho=r sin(\theta)','FontSize',22);
xlabel('z=r cos(\theta)','FontSize',22);



end
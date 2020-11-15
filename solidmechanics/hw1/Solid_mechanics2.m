clc;clear;
%Resize the grid;detla.larger delta;larger divided.
detla = 1;
%unit (km)
x=-75:1.0/detla:250;y=0:1.0/detla:300;
[X,Y]=meshgrid(x,y);[m,n]=size(X);
%standard unit
X=1000*X;Y=1000*Y;
%We set the U miu and theta
U_true=5;U=U_true*10^(-2)/(3600*(24*365+6));miu=10^20;theta=pi/6;
%continental crust
%(u=v=0,y=0,x>0;u=3^(1/2)/2*U,v=1/2*U,x=3^(1/2)*y)
Cc=[0 -1 -1 0;
    1 0 0 0;
    0 -1 -3/4 -(pi/6+3^(1/2)/4);
    1 0 (pi/6-3^(1/2)/4) -1/4];
C_uv=[0 0 3^(1/2)/2*U 1/2*U]';
Para_c=Cc\C_uv;
Aa=Para_c(1);Ba=Para_c(2);Ca=Para_c(3);Da=Para_c(4);
%Oceanic crust
%(v=0,u=U,y=0,x<0;u=3^(1/2)/2*U,v=1/2*U,x=3^(1/2)*y)
Oc=[0 -1 -1 -pi;
     1 0 pi 0;
    0 -1 -3/4 -(pi/6+3^(1/2)/4);
    1 0 (pi/6-3^(1/2)/4) -1/4];
O_uv=[U 0 3^(1/2)/2*U 1/2*U]';
Para_o=Oc\O_uv;
Ao=Para_o(1);Bo=Para_o(2);Co=Para_o(3);Do=Para_o(4);

u=zeros(m,n);v=u;P=u;
% The distribution of P;U;V in the whole space
for i=1:m
    for j=1:n
        if X(i,j)>(Y(i,j)*3^(1/2))% arc corner
            u(i,j)=-Ba-Da*atan(Y(i,j)/X(i,j))+(Ca*X(i,j)+Da*Y(i,j))*(-X(i,j)/(X(i,j)^2+Y(i,j)^2));
            v(i,j)=Aa+Ca*atan(Y(i,j)/X(i,j))+(Ca*X(i,j)+Da*Y(i,j))*(-Y(i,j)/(X(i,j)^2+Y(i,j)^2));
            P(i,j)=-2*miu*(Ca*X(i,j)+Da*Y(i,j))/(X(i,j)^2+Y(i,j)^2);
        else % oceanic corner
            if X(i,j)>0
                u(i,j)=-Bo-Do*atan(Y(i,j)/X(i,j))+(Co*X(i,j)+Do*Y(i,j))*(-X(i,j)/(X(i,j)^2+Y(i,j)^2));
                v(i,j)=Ao+Co*atan(Y(i,j)/X(i,j))+(Co*X(i,j)+Do*Y(i,j))*(-Y(i,j)/(X(i,j)^2+Y(i,j)^2));
            else % the matlb function has some problem; if dont't add pi, there will be no consistant.
                u(i,j)=-Bo-Do*(atan(Y(i,j)/X(i,j))+pi)+(Co*X(i,j)+Do*Y(i,j))*(-X(i,j)/(X(i,j)^2+Y(i,j)^2));
                v(i,j)=Ao+Co*(atan(Y(i,j)/X(i,j))+pi)+(Co*X(i,j)+Do*Y(i,j))*(-Y(i,j)/(X(i,j)^2+Y(i,j)^2));
            end
            P(i,j)=-2*miu*(Co*X(i,j)+Do*Y(i,j))/(X(i,j)^2+Y(i,j)^2);
        end
    end
end
%%
tao_xx=zeros(m,n);tao_yy=tao_xx;tao_xy=tao_xx;tao_yx=tao_xx;
max_shear= tao_xx;sigma1=tao_xx;sigma2=tao_xx;
S_ij=zeros(2,2);%The deviatoric stress tensor
%if here are NAN; turn to 0.
P(isnan(P))=0;
u(isnan(u))=0;
v(isnan(v))=0;
%Take the partial of ux and u sub y
[ux,uy]=gradient(u,1000.,1000.);[vx,vy]=gradient(v,1000.,1000.);
%calculate the tao
tao_xx=2*miu*ux*10^(-6);tao_yy=2*miu*vy*10^(-6);
tao_xy=miu*(uy+vx)*10^(-6);tao_yx=miu*(uy+vx)*10^(-6);
P=P*10^(-6);
%Turn the velocity back to cm/year
u=u*(10^2*3600*(24*365+6));
v=v*(10^2*3600*(24*365+6));

%Principal stress and shear stress
for i=1:m
    for j=1:n
        S_ij=[(P(i,j)-tao_xx(i,j)) tao_xy(i,j);tao_yx(i,j) (P(i,j)-tao_yy(i,j))];
        sigma=eig(S_ij);
        sigma1(i,j)=sigma(1);
        sigma2(i,j)=sigma(2);
        max_shear(i,j)=abs(sigma1(i,j)-sigma2(i,j))/2;
        %the direction of large one
        theta(i,j)=1./2*atan(2*tao_xy(i,j)/(tao_xx(i,j)-tao_yy(i,j)));
        sigma1_x1(i,j)=-sigma1(i,j)/2*cos(theta(i,j));
        sigma1_y1(i,j)=-sigma1(i,j)/2*sin(theta(i,j));
        sigma1_x2(i,j)=+sigma1(i,j)/2*cos(theta(i,j));
        sigma1_y2(i,j)=+sigma1(i,j)/2*sin(theta(i,j));
        %the direction of smaller one
        theta(i,j)=theta(i,j)+pi/2;
        sigma2_x1(i,j)=-sigma2(i,j)/2*cos(theta(i,j));
        sigma2_y1(i,j)=-sigma2(i,j)/2*sin(theta(i,j));
        sigma2_x2(i,j)=+sigma2(i,j)/2*cos(theta(i,j));
        sigma2_y2(i,j)=+sigma2(i,j)/2*sin(theta(i,j));
    end
end
%% 

%figure 
figure(1);
pcolor(X/1000,Y/1000,P);
shading interp;
axis ij;
colormap(figure(1),hot);colorbar;
h=colorbar;
set(gca,'CLim',[-50,20]);
hold on;
set(get(h,'Title'),'string','MPa');
xlabel('km');
ylabel('km');
hold on;
title('Velocity and pressure fields wth U=5 cm/year  and theta=30бу');
step=20;
scale=4;
dx=-75:step:250;
dy=0:step:300;
[Dx,Dy]=meshgrid(dx,dy);
fenbu_u=u(1:detla*step:end,1:detla*step:end);
fenbu_v=v(1:detla*step:end,1:detla*step:end);
quiver(Dx,Dy,fenbu_u*scale,fenbu_v*scale,'autoscale','off','color','[0.5 0.5 0.5]'); 
hold on;
quiver(40,12,U_true*scale*2,0*scale*2,'autoscale','off','color','w');
text(50,5,'5 cm/year','color','w')
hold off;  
%%
figure(2);
pcolor(X/1000,Y/1000,max_shear);
shading interp;axis ij;
colormap(figure(2),parula);colorbar;
h=colorbar;
set(gca,'CLim',[0,30]);
hold on;
set(get(h,'Title'),'string','MPa');
xlabel('km');ylabel('km');
hold on;
title({'Principal stress and maximum shear stress'});
% cross
step=15;
scale=4; 
dx=-25:step:300;
dy=0:step:300;
[DX,DY]=meshgrid(dx,dy);
Cross1=sigma1(1:detla*step:end,1:detla*step:end);
Cross2=sigma2(1:detla*step:end,1:detla*step:end);
[m,n]=size(Cross1);
%max
Cross1_x1=sigma1_x1(1:detla*step:end,1:detla*step:end);
Cross1_x2=sigma1_x2(1:detla*step:end,1:detla*step:end);
Cross1_y1=sigma1_y1(1:detla*step:end,1:detla*step:end);
Cross1_y2=sigma1_y2(1:detla*step:end,1:detla*step:end);
%min
Cross2_x1=sigma2_x1(1:detla*step:end,1:detla*step:end);
Cross2_x2=sigma2_x2(1:detla*step:end,1:detla*step:end);
Cross2_y1=sigma2_y1(1:detla*step:end,1:detla*step:end);
Cross2_y2=sigma2_y2(1:detla*step:end,1:detla*step:end);
% 
for i=1:m
    for j=1:n
        x1=[DX(i,j)+Cross1_x1(i,j)*scale DX(i,j)+Cross1_x2(i,j)*scale];
        y1=[DY(i,j)+Cross1_y1(i,j)*scale DY(i,j)+Cross1_y2(i,j)*scale];
        x2=[DX(i,j)+Cross2_x1(i,j)*scale DX(i,j)+Cross2_x2(i,j)*scale];
        y2=[DY(i,j)+Cross2_y1(i,j)*scale DY(i,j)+Cross2_y2(i,j)*scale];
        if abs(x1(1)-x1(2))>80||abs(x2(1)-x2(2))>80
            continue
        end
           
        if Cross1(i,j)>=0
            plot(x1,y1,'r-');
        else 
            plot(x1,y1,'w-');
        end
        hold on;
        
      
        if Cross2(i,j)>=0
            plot(x2,y2,'r-');
        else 
            plot(x2,y2,'w-');
        end
        hold on;
    end
end
hold on;
% scale
x=[230, 230+scale*2];y=[252,252];
plot(x,y,'y-');
text(230,245,'2.5Mpa','color','w');
hold off;

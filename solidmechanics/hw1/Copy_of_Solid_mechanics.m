clc;clear;

%构建�??标系和网格（�?��?/km）
x=-25:1:300;
y=0:1:300; 
[X,Y]=meshgrid(x,y);
[m,n]=size(X);
%输入速度�?�数并转化�?��?（m/s）
U_true=5;%(cm/year)
U=U_true*10^(-2)/(3600*24*365);
%输入粘度拉梅常数（Pa*s）
miu=10^20;
%输入角度
theta1=pi/6;

%开始进行求解速度和压强分布
%u和v的公�?由�?�考书�?获�?�
%对于arc corner则有下�?�求解过程
   %将u和v的函数表达�?�?
   %(u=v=0,y=0,x>0;u=3^(1/2)/2*U,v=1/2*U,x=3^(1/2)*y)
xs=[0 -1 -1 0;1 0 0 0;
    0 -1 -3/2 -(pi/6+3^(1/2)/4);1 0 (pi/6-3^(1/2)/4) -1/4];%系数矩阵
rt=[0 0 3^(1/2)/2*U 1/2*U]';%结果矩阵
cs=xs\rt;
Aa=cs(1);Ba=cs(2);Ca=cs(3);Da=cs(4);
%�?��?��?�以求得ocean corner相关分布
   %(v=0,u=U,y=0,x<0;u=3^(1/2)/2*U,v=1/2*U,x=3^(1/2)*y)
xs1=[0 -1 -1 -pi;1 0 pi 0;
    0 -1 -3/2 -(pi/6+3^(1/2)/4);1 0 (pi/6-3^(1/2)/4) -1/4];
rt1=[U 0 3^(1/2)/2*U 1/2*U]';
cs1=xs1\rt1;
Ao=cs1(1);Bo=cs1(2);Co=cs1(3);Do=cs1(4);


X=1000*X;Y=1000*Y;

u=zeros(m,n);v=zeros(m,n);
P=zeros(m,n);

for i=1:m
    for j=1:n
        if X(i,j)>(Y(i,j)*3^(1/2))%在arc corner一边
            u(i,j)=-Ba-Da*atan(Y(i,j)/X(i,j))+(Ca*X(i,j)+Da*Y(i,j))*(-X(i,j)/(X(i,j)^2+Y(i,j)^2));
            v(i,j)=Aa+Ca*atan(Y(i,j)/X(i,j))+(Ca*X(i,j)+Da*Y(i,j))*(-Y(i,j)/(X(i,j)^2+Y(i,j)^2));
            P(i,j)=-2*miu*(Ca*X(i,j)+Da*Y(i,j))/(X(i,j)^2+Y(i,j)^2);
        else %在oceanic corner 一边
            if X(i,j)>0
                u(i,j)=-Bo-Do*atan(Y(i,j)/X(i,j))+(Co*X(i,j)+Do*Y(i,j))*(-X(i,j)/(X(i,j)^2+Y(i,j)^2));
                v(i,j)=Ao+Co*atan(Y(i,j)/X(i,j))+(Co*X(i,j)+Do*Y(i,j))*(-Y(i,j)/(X(i,j)^2+Y(i,j)^2));
            else
                u(i,j)=-Bo-Do*(atan(Y(i,j)/X(i,j))+pi)+(Co*X(i,j)+Do*Y(i,j))*(-X(i,j)/(X(i,j)^2+Y(i,j)^2));
                v(i,j)=Ao+Co*(atan(Y(i,j)/X(i,j))+pi)+(Co*X(i,j)+Do*Y(i,j))*(-Y(i,j)/(X(i,j)^2+Y(i,j)^2));
            end
            P(i,j)=-2*miu*(Co*X(i,j)+Do*Y(i,j))/(X(i,j)^2+Y(i,j)^2);
        end
    end
end



tao_xx=zeros(m,n);
tao_yy=zeros(m,n);
tao_xy=zeros(m,n);
tao_yx=zeros(m,n);
max_shear=zeros(m,n);%最大剪应力
sigma1=zeros(m,n);
sigma2=zeros(m,n);%两个主应力

A=zeros(2,2);

P(isnan(P))=0;
u(isnan(u))=0;
v(isnan(v))=0;

[ux,uy]=gradient(u,1000.,1000.);
[vx,vy]=gradient(v,1000.,1000.);

tao_xx=2*miu*ux;
tao_yy=2*miu*vy;
tao_xy=miu*(uy+vx);
tao_yx=miu*(uy+vx);

tao_xx=tao_xx*10^(-6);
tao_yy=tao_yy*10^(-6);
tao_xy=tao_xy*10^(-6);
tao_yx=tao_yx*10^(-6);
P=P*10^(-6);
%将分�?速度�?��?�?新转化为cm/year
u=u*10^2*365*24*3600;
v=v*10^2*365*24*3600;

%根�?��??应力张�?求�?�主应力和最大剪应力
for i=1:m
    for j=1:n
        A=[tao_xx(i,j) tao_xy(i,j);tao_yx(i,j) tao_yy(i,j)];
        sigma=eig(A);
        sigma1(i,j)=sigma(1);
        sigma2(i,j)=sigma(2);
        max_shear(i,j)=abs(sigma1(i,j)-sigma2(i,j))/2;%（值�?�能有问题）
        %算出最大主应力的方�?�
        theta(i,j)=1./2*atan(2*tao_xy(i,j)/(tao_xx(i,j)-tao_yy(i,j)));
        sigma1_x1(i,j)=-sigma1(i,j)/2*cos(theta(i,j));
        sigma1_y1(i,j)=-sigma1(i,j)/2*sin(theta(i,j));
        sigma1_x2(i,j)=+sigma1(i,j)/2*cos(theta(i,j));
        sigma1_y2(i,j)=+sigma1(i,j)/2*sin(theta(i,j));
        %算出最�?主应力的方�?�
        theta(i,j)=theta(i,j)+pi/2;
        sigma2_x1(i,j)=-sigma2(i,j)/2*cos(theta(i,j));
        sigma2_y1(i,j)=-sigma2(i,j)/2*sin(theta(i,j));
        sigma2_x2(i,j)=+sigma2(i,j)/2*cos(theta(i,j));
        sigma2_y2(i,j)=+sigma2(i,j)/2*sin(theta(i,j));
    end
end
%% 

%开始画图
figure(1);
pcolor(X/1000,Y/1000,P);%画出压强
shading interp;%（使颜色平滑分布）
axis ij;%（使y值从上到下增加,�??标轴改为矩阵模�?）
h=colorbar;
set(gca,'CLim',[-50,20]);%(调整colorbar)
hold on;
set(get(h,'Title'),'string','MPa');%(改colorbar�??字)
xlabel('km');
ylabel('km');
hold on;
title('速度和压强分布图 | U=5 cm/year theta=30度');
step=20;
scale=2;
dx=-25:step:300;
dy=0:step:300;
[Dx,Dy]=meshgrid(dx,dy);
fenbu_u=u(1:step:end,1:step:end);
fenbu_v=v(1:step:end,1:step:end);%画出速度分布图
quiver(Dx,Dy,fenbu_u*scale,fenbu_v*scale,'autoscale','off','color','k'); 
%（画出矢�?分布箭头，采用给定的比例尺，所以关闭自动分�?的比例)
hold on;
quiver(260,12,U_true*scale*2,0*scale*2,'autoscale','off','color','k');
%(给定比例尺大�?)
text(265,5,'5 cm/year','color','k')%(定义字符串�?�在�?个�?置添加文字)
hold off;  

%画出主应力和最大剪应力方�?�
figure(2);
%画出最大剪应力
pcolor(X/1000,Y/1000,max_shear);
shading interp;axis ij;
h=colorbar;
set(gca,'CLim',[0,30]);
hold on;
set(get(h,'Title'),'string','MPa');
xlabel('km');ylabel('km');
hold on;
title({'主应力和最大剪应力'});
% 主应力�??字
step=15;
scale=3; 
dx=-25:step:300;
dy=0:step:300;
[DX,DY]=meshgrid(dx,dy);
Cross1=sigma1(1:step:end,1:step:end);
Cross2=sigma2(1:step:end,1:step:end);
[m,n]=size(Cross1);
% 最大主应力
Cross1_x1=sigma1_x1(1:step:end,1:step:end);
Cross1_x2=sigma1_x2(1:step:end,1:step:end);
Cross1_y1=sigma1_y1(1:step:end,1:step:end);
Cross1_y2=sigma1_y2(1:step:end,1:step:end);
% 最�?主应力 
Cross2_x1=sigma2_x1(1:step:end,1:step:end);
Cross2_x2=sigma2_x2(1:step:end,1:step:end);
Cross2_y1=sigma2_y1(1:step:end,1:step:end);
Cross2_y2=sigma2_y2(1:step:end,1:step:end);
% 画出主应力�??字
for i=1:m
    for j=1:n
        x1=[DX(i,j)+Cross1_x1(i,j)*scale DX(i,j)+Cross1_x2(i,j)*scale];
        y1=[DY(i,j)+Cross1_y1(i,j)*scale DY(i,j)+Cross1_y2(i,j)*scale];
        x2=[DX(i,j)+Cross2_x1(i,j)*scale DX(i,j)+Cross2_x2(i,j)*scale];
        y2=[DY(i,j)+Cross2_y1(i,j)*scale DY(i,j)+Cross2_y2(i,j)*scale];
        if abs(x1(1)-x1(2))>80||abs(x2(1)-x2(2))>80
            continue
        end
        % 画出最大主应力
        if Cross1(i,j)>=0
            plot(x1,y1,'r-');
        else 
            plot(x1,y1,'w-');
        end
        hold on;
        
        % 画出最�?主应力

        if Cross2(i,j)>=0
            plot(x2,y2,'r-');
        else 
            plot(x2,y2,'w-');
        end
        hold on;
    end
end
hold on;
% 画出比例尺
x=[275, 275+scale*2];y=[12,12];
plot(x,y,'y-');
text(275,5,'2.5Mpa','color','y');
hold off;

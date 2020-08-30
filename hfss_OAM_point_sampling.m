%%此程序为从仿真软件中获得得到的幅度相位数据文件采样模态谱图绘制
%仅供参考学习 Ling-Jun Yang    2020.8.25
% clc ,clear
clc , close all
E1 = importdata('single_squarloop9GHzmode1.fld');  %d导入矢量电场cvc类型
% E2 = importdata('single_squarloop9GHzinc.fld');


%以下的部分代码是用来将数据文件转化为矩阵
E=E1.data;
% Ein=E2.data;
Ein=zeros(size(E));
[MM,NN]=size(E); %读出导入矩阵的大小（九行mm*nn列）

mm =201;%行维数，必须保证为整数，数据矩阵 与实际的数据文件对应
nn =201;%列维数,     mm行，nn列
oo=8;%页
efftotal=ones(1,oo);
efftotal1=ones(1,oo);

    
Ex=zeros(mm,nn,oo);
Ey =zeros(mm,nn,oo);

y=1;%一个设置的变量，便于循环的数据点选取
for p = 1:mm
    for q=1:nn
 lie1=(E(y:y+oo-1,4)-Ein(y:y+oo-1,4))+1j.*(E(y:y+oo-1,5)-Ein(y:y+oo-1,5));  %用于导入相位数据
 Ex(p,q,:)=lie1;
 
 lie2=E(y:y+oo-1,6)-Ein(y:y+oo-1,6)+1j.*E(y:y+oo-1,7)-1j.*Ein(y:y+oo-1,7); %用于导入幅度数据
 Ey(p,q,:)=lie2;
 y=y+oo; %迭代，产生下一个数据的起始点
    end
end  
%提取圆极化：
El=zeros(mm,nn,oo);
Er =zeros(mm,nn,oo);
t=1/sqrt(2)*[1,-1j;1,1j];%t为转换矩阵
for p = 1:mm
    for q = 1:nn
        for o=1:oo
      El(p,q,o)  = 1/sqrt(2)*( Ex(p,q,o)-1j*Ey(p,q,o)) ; %线极化变圆计划
      Er(p,q,o)  = 1/sqrt(2)*( Ex(p,q,o)+1j*Ey(p,q,o));
        end
    end
end  



% Ex=Ex-ones(size(Ex));
%获得矩阵形式的幅度相位分布
phase_x=angle(Ex)/pi*180;
phase_l=angle(El)/pi*180;
phase_r=angle(Er)/pi*180;
% phase_x=atan(real(Ey)./real(Ex));
mag_E=sqrt(abs(Ex).^2+abs(Ey).^2);
mag_El=abs(El);
mag_Er=abs(Er);
 AX=linspace(-0.4,0.4,201);
AY=linspace(-0.4,0.4,201);



n_mag_A1=ones(oo,15);
n_mag_B1=ones(oo,15);
for o=1:oo
maxvalue=max(max(mag_El(:,:,o)));    %求出幅度最大值
centre=101; %由剖分的网格步长和网格大小以及多少确定的网格中心点，相对位置centre
[row,col]=find(maxvalue==mag_El(:,:,o));%求出幅度最大值位置
s_R=sqrt((row-centre)^2+(col-centre)^2);
if s_R>100
    s_R=100;
end
s_N=100;
dtheta=2*pi/s_N;
s_theta=dtheta:dtheta:2*pi;%角度分布数组,注意此处,表示从dtheta度开始取样
s_x =s_R.* cos(s_theta); %阵元坐标x轴  数组
s_y =s_R.* sin (s_theta);%阵元坐标y轴  数组
% %以下是用来计算坐标的
s_x=s_x+centre; %301是原始MATLAB代码中场观察截面的中心点
s_y=s_y+centre;
s_x=round(s_x);  %四舍五入取整数
s_y=round(s_y);
N_mode=15;
A=ones(1,N_mode);%左旋分量模态谱
B=ones(1,N_mode);%右旋分量模态谱
for k=1:N_mode
    mode=0+0*1j;
    mode1=0+0*1j;
    for i=1:s_N
        mode=mode+dtheta*El(s_x(i),s_y(i),o)*exp(-1j*(k-8)*s_theta(i));
        mode1=mode1+dtheta*Er(s_x(i),s_y(i),o)*exp(-1j*(k-8)*s_theta(i));
    end
    A(k)=mode;
    B(k)=mode1;
end
mag_A=abs(A);
mag_B=abs(B);
n_mag_A=mag_A/max(mag_A);
n_mag_B=mag_B/max(mag_A);
eff=1/(sum(n_mag_A)+sum(n_mag_B));
eff1=1/(n_mag_B(8)+1);
efftotal(o)=eff;
efftotal1(o)=eff1;
n_mag_A1(o,:)=mag_A/max(mag_A);
n_mag_B1(o,:)=mag_B/max(mag_A);

end
a=max(efftotal1);
% o=find(efftotal1==a)
o=6;
%%%% 绘图

 figure(1)%OAM谱
 subplot(1,2,1)
 stem(-7:7,n_mag_A1(o,:));
% subplot(2,1,1);stem(0:length(X1)-1,magX);
 xlabel('l');ylabel('(2)LH');title('magnitude part')
  subplot(1,2,2)
 stem(-7:7,n_mag_B1(o,:));
  ylim([0 1])
% subplot(2,1,1);stem(0:length(X1)-1,magX);
 xlabel('l');ylabel('(2)RH');title('magnitude part')

maxcbar=2;

%%%%%%画左旋右旋
figure(3)
imagesc(AX,AY,mag_El(:,:,o))
colormap(jet)
axis xy; %此处是将坐标轴调整正
colorbar
xlabel('x-axis (m)','Fontname', 'Times New Roman','FontSize',8)
ylabel('y-axis (m)','Fontname', 'Times New Roman','FontSize',8)
zlabel('z-axis (m)','Fontname', 'Times New Roman','FontSize',8)
set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1,'FontWeight','bold')
title('E-field intensity')
h=colorbar
set(h,'FontName','Times New Roman', 'FontSize',16, 'FontWeight', 'bold') %colorbar 字体设置 加粗
caxis([0,maxcbar]);
hold on 


 figure(4) 
% subplot(1,2,2),
imagesc(AX,AY,phase_l(:,:,o))
% colormap(jet)
 colormap(hsv)
axis xy; %此处是将坐标轴调整正
colorbar;
xlabel('x-axis (m)','Fontname', 'Times New Roman','FontSize',8)
ylabel('y-axis (m)','Fontname', 'Times New Roman','FontSize',8)
zlabel('z-axis (m)','Fontname', 'Times New Roman','FontSize',8)
set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1,'FontWeight','bold')
title('E-field phase')
hh=colorbar;
set(hh,'FontName','Times New Roman', 'FontSize',16, 'FontWeight', 'bold') %colorbar 字体设置 加粗
set(get(hh,'Title'),'string','Phase (deg)');
hold on 

figure(5)
imagesc(AX,AY,mag_Er(:,:,o))
colorbar;
colormap(jet)
axis xy; %此处是将坐标轴调整正
xlabel('x-axis (m)','Fontname', 'Times New Roman','FontSize',8)
ylabel('y-axis (m)','Fontname', 'Times New Roman','FontSize',8)
zlabel('z-axis (m)','Fontname', 'Times New Roman','FontSize',8)
set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1,'FontWeight','bold')
title('E-field intensity')
h=colorbar
set(h,'FontName','Times New Roman', 'FontSize',16, 'FontWeight', 'bold') %colorbar 字体设置 加粗
caxis([0,maxcbar]);
hold on 
% view(0,-90);

figure(6) 
imagesc(AX,AY,phase_r(:,:,o))
colormap(hsv)
colorbar;
xlabel('x-axis (m)','Fontname', 'Times New Roman','FontSize',8)
ylabel('y-axis (m)','Fontname', 'Times New Roman','FontSize',8)
zlabel('z-axis (m)','Fontname', 'Times New Roman','FontSize',8)
set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1,'FontWeight','bold')
title('E-field phase')
hh=colorbar;
set(hh,'FontName','Times New Roman', 'FontSize',16, 'FontWeight', 'bold') %colorbar 字体设置 加粗
set(get(hh,'Title'),'string','Phase (deg)');
hold on 






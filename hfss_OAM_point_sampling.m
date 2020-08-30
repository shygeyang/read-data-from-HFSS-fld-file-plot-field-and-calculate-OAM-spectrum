%%�˳���Ϊ�ӷ�������л�õõ��ķ�����λ�����ļ�����ģ̬��ͼ����
%�����ο�ѧϰ Ling-Jun Yang    2020.8.25
% clc ,clear
clc , close all
E1 = importdata('single_squarloop9GHzmode1.fld');  %d����ʸ���糡cvc����
% E2 = importdata('single_squarloop9GHzinc.fld');


%���µĲ��ִ����������������ļ�ת��Ϊ����
E=E1.data;
% Ein=E2.data;
Ein=zeros(size(E));
[MM,NN]=size(E); %�����������Ĵ�С������mm*nn�У�

mm =201;%��ά�������뱣֤Ϊ���������ݾ��� ��ʵ�ʵ������ļ���Ӧ
nn =201;%��ά��,     mm�У�nn��
oo=8;%ҳ
efftotal=ones(1,oo);
efftotal1=ones(1,oo);

    
Ex=zeros(mm,nn,oo);
Ey =zeros(mm,nn,oo);

y=1;%һ�����õı���������ѭ�������ݵ�ѡȡ
for p = 1:mm
    for q=1:nn
 lie1=(E(y:y+oo-1,4)-Ein(y:y+oo-1,4))+1j.*(E(y:y+oo-1,5)-Ein(y:y+oo-1,5));  %���ڵ�����λ����
 Ex(p,q,:)=lie1;
 
 lie2=E(y:y+oo-1,6)-Ein(y:y+oo-1,6)+1j.*E(y:y+oo-1,7)-1j.*Ein(y:y+oo-1,7); %���ڵ����������
 Ey(p,q,:)=lie2;
 y=y+oo; %������������һ�����ݵ���ʼ��
    end
end  
%��ȡԲ������
El=zeros(mm,nn,oo);
Er =zeros(mm,nn,oo);
t=1/sqrt(2)*[1,-1j;1,1j];%tΪת������
for p = 1:mm
    for q = 1:nn
        for o=1:oo
      El(p,q,o)  = 1/sqrt(2)*( Ex(p,q,o)-1j*Ey(p,q,o)) ; %�߼�����Բ�ƻ�
      Er(p,q,o)  = 1/sqrt(2)*( Ex(p,q,o)+1j*Ey(p,q,o));
        end
    end
end  



% Ex=Ex-ones(size(Ex));
%��þ�����ʽ�ķ�����λ�ֲ�
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
maxvalue=max(max(mag_El(:,:,o)));    %����������ֵ
centre=101; %���ʷֵ����񲽳��������С�Լ�����ȷ�����������ĵ㣬���λ��centre
[row,col]=find(maxvalue==mag_El(:,:,o));%����������ֵλ��
s_R=sqrt((row-centre)^2+(col-centre)^2);
if s_R>100
    s_R=100;
end
s_N=100;
dtheta=2*pi/s_N;
s_theta=dtheta:dtheta:2*pi;%�Ƕȷֲ�����,ע��˴�,��ʾ��dtheta�ȿ�ʼȡ��
s_x =s_R.* cos(s_theta); %��Ԫ����x��  ����
s_y =s_R.* sin (s_theta);%��Ԫ����y��  ����
% %�������������������
s_x=s_x+centre; %301��ԭʼMATLAB�����г��۲��������ĵ�
s_y=s_y+centre;
s_x=round(s_x);  %��������ȡ����
s_y=round(s_y);
N_mode=15;
A=ones(1,N_mode);%��������ģ̬��
B=ones(1,N_mode);%��������ģ̬��
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
%%%% ��ͼ

 figure(1)%OAM��
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

%%%%%%����������
figure(3)
imagesc(AX,AY,mag_El(:,:,o))
colormap(jet)
axis xy; %�˴��ǽ������������
colorbar
xlabel('x-axis (m)','Fontname', 'Times New Roman','FontSize',8)
ylabel('y-axis (m)','Fontname', 'Times New Roman','FontSize',8)
zlabel('z-axis (m)','Fontname', 'Times New Roman','FontSize',8)
set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1,'FontWeight','bold')
title('E-field intensity')
h=colorbar
set(h,'FontName','Times New Roman', 'FontSize',16, 'FontWeight', 'bold') %colorbar �������� �Ӵ�
caxis([0,maxcbar]);
hold on 


 figure(4) 
% subplot(1,2,2),
imagesc(AX,AY,phase_l(:,:,o))
% colormap(jet)
 colormap(hsv)
axis xy; %�˴��ǽ������������
colorbar;
xlabel('x-axis (m)','Fontname', 'Times New Roman','FontSize',8)
ylabel('y-axis (m)','Fontname', 'Times New Roman','FontSize',8)
zlabel('z-axis (m)','Fontname', 'Times New Roman','FontSize',8)
set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1,'FontWeight','bold')
title('E-field phase')
hh=colorbar;
set(hh,'FontName','Times New Roman', 'FontSize',16, 'FontWeight', 'bold') %colorbar �������� �Ӵ�
set(get(hh,'Title'),'string','Phase (deg)');
hold on 

figure(5)
imagesc(AX,AY,mag_Er(:,:,o))
colorbar;
colormap(jet)
axis xy; %�˴��ǽ������������
xlabel('x-axis (m)','Fontname', 'Times New Roman','FontSize',8)
ylabel('y-axis (m)','Fontname', 'Times New Roman','FontSize',8)
zlabel('z-axis (m)','Fontname', 'Times New Roman','FontSize',8)
set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1,'FontWeight','bold')
title('E-field intensity')
h=colorbar
set(h,'FontName','Times New Roman', 'FontSize',16, 'FontWeight', 'bold') %colorbar �������� �Ӵ�
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
set(hh,'FontName','Times New Roman', 'FontSize',16, 'FontWeight', 'bold') %colorbar �������� �Ӵ�
set(get(hh,'Title'),'string','Phase (deg)');
hold on 






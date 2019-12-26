clc
clear;
%%      ��������2��ʱ���֣��ռ�2*6 �ײ�֡�
%%%%%%%%%%%%%%%%%%%%%%�ο����£�����ȫƥ������ձ߽�����Ӧ���о�������ά��ȣ�2013
%%%%%%%��������
Nx = 500; Nz = 600; Nt = 2300;
dt =0.001; dx = 10;dz =10;
v1 = 3000;  %�ٶ�
v2 = 2500;
f = 20;    % ��Դ������Ƶ
% a = [1.5 -0.15 0.011111];   %���ϵ��
a = [1.66667 -0.238095 0.0396825 -0.00496032 0.000317460];   %���ϵ��
%%
PML_n = 100;  %���ձ߽����
n=5;  % �ռ�n*2�ײ��
sum_Nx = PML_n*2+n*2+Nx;
sum_Nz = PML_n*2+n*2+Nz;
v=zeros(sum_Nx,sum_Nz)+v1;
v(1:sum_Nz/2,:) = v2;
s_x =sum_Nx/2+200;s_z =sum_Nz/2;

%%
A1 =zeros(sum_Nx,sum_Nz);       %˥������
A2 =zeros(sum_Nx,sum_Nz);       %˥������
A3 =zeros(sum_Nx,sum_Nz);       %˥������
A4 =zeros(sum_Nx,sum_Nz);       %˥������
A =zeros(sum_Nx,sum_Nz);       %˥������
u = zeros(sum_Nx,sum_Nz,Nt);   % ����ֵ

traces = 21;
start_traces = PML_n;
traces_interval = 30;
seis_wave = zeros(traces,Nt);
%%%%%%%%%%%%%%%%%%%%%%%%%% ���ٶ�ģ��
figure(2)
imagesc(v(n+1:sum_Nx-n,n+1:sum_Nz-n))
hold on;
scatter([start_traces:traces_interval:start_traces+traces_interval*(traces-1)],(1+n+PML_n)*ones(1,traces),200,'rv','filled')
hco = colorbar;
set(get(hco,'ylabel'),'string','velocity/m*s^-^1')
set(gca,'xAxislocation','top');
xlabel('x','fontname','times new roman','fontsize',30); ylabel('z','fontname','times new roman','fontsize',30);
set(gca,'fontname','times new roman','fontsize',30)
set(gcf,'unit','centimeters','position',[5 3 33 20])
saveas (gcf,['two_layers_velocity','.png']);
% if v*(dt/dx)<=1
%     disp('�����ȶ�������');
% else
%     disp('�������ȶ�������');
% end
% for k = 2: 100
%     u(k) = (1-2*(pi*f*((k-50)*dt)).^2).*exp(-(pi*f*((k-50)*dt)).^2);
% end
% plot(u);
%%        ��ʼ���㲨��,��ѡ����ȫƥ������ձ߽磨PML��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ָ��˥������1
% R_coeff = 0.0001;
% for i=n+1:sum_Nx-n
%     for j=n+1:sum_Nz-n
%         if j>n&&j<=PML_n+n
%             %                 A1(i,j) = B*(1-cos(pi*(PML_n+n+1-j)/(2*PML_n)));  %������
%             A1(i,j) = 3*v(i,j)/(2*PML_n)*log(1/R_coeff)*((PML_n+n+1-j)/PML_n)^2;  %������
%         elseif j>n+PML_n+Nx&&j<=n+PML_n+Nx+PML_n
%             %                 A2(i,j) = B*(1-cos(pi*(j-(n+PML_n+Nz))/(2*PML_n)));  %������
%             A2(i,j) = 3*v(i,j)/(2*PML_n)*log(1/R_coeff)*((j-(n+PML_n+Nz))/PML_n)^2;  %������
%         end
%         
%         if i>n&&i<=PML_n+n
%             %                 A3(i,j) = B*(1-cos(pi*(PML_n+n+1-i)/(2*PML_n)));  %������
%             A3(i,j) = 3*v(i,j)/(2*PML_n)*log(1/R_coeff)*((PML_n+n+1-i)/PML_n)^2;  %������
%         elseif i>n+PML_n+Nx&&i<=n+PML_n+Nx+PML_n
%             %                 A4(i,j) = B*(1-cos(pi*(i-(n+PML_n+Nx))/(2*PML_n)));  %������
%             A4(i,j) = 3*v(i,j)/(2*PML_n)*log(1/R_coeff)*((i-(n+PML_n+Nx))/PML_n)^2;  %������
%         end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%������˥������2
B = 100;  % ˥����������
for i=n+1:sum_Nx-n
    for j=n+1:sum_Nz-n
        if j>n&&j<=PML_n+n
            A1(i,j) = B*(1-cos(pi*(PML_n+n+1-j)/(2*PML_n)));  %������
        elseif j>n+PML_n+Nz&&j<=n+PML_n+Nz+PML_n
            A2(i,j) = B*(1-cos(pi*(j-(n+PML_n+Nz))/(2*PML_n)));  %������
        end
        
        if i>n&&i<=PML_n+n
            A3(i,j) = B*(1-cos(pi*(PML_n+n+1-i)/(2*PML_n)));  %������
        elseif i>n+PML_n+Nx&&i<=n+PML_n+Nx+PML_n
            A4(i,j) = B*(1-cos(pi*(i-(n+PML_n+Nx))/(2*PML_n)));  %������
        end
    end
end
A=A1+A2+A3+A4;
clear A1 A2 A3 A4
%%%%%%%%%%%%%%%%%%%%%%%%  ��˥�����ӷֲ�ͼ
figure(3)
imagesc(A(n+1:sum_Nx-n,n+1:sum_Nz-n))
shading interp;
% hco = colorbar;
% set(get(hco,'ylabel'),'string','velocity/m*s^-^1')
set(gca,'Ydir','reverse','xAxislocation','top');
xlabel('x','fontname','times new roman','fontsize',30); ylabel('z','fontname','times new roman','fontsize',30);
set(gca,'fontname','times new roman','fontsize',30)
set(gcf,'unit','centimeters','position',[5 3 33 20])
saveas (gcf,['˥�����ӷֲ�ͼ','.png']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���㲨��ֵ
for k = 2 :Nt
    %     u(s_i+PML_n+n,s_j+PML_n+n,k) = (1-2*(pi*f*((k-30)*dt)).^2).*exp(-(pi*f*((k-30)*dt)).^2);  %�� ��201��201���������Դʱ�亯��
    for i = n+1: sum_Nx-n
        for j = n+1: sum_Nz-n
            if i==s_x && j ==s_z
                u(i,j,k) = (1-2*(pi*f*((k-50)*dt)).^2).*exp(-(pi*f*((k-50)*dt)).^2);  %�� ��201��201���������Դʱ�亯��
            else
                coe1 = (2-A(i,j)*A(i,j)*dt*dt)/(1+A(i,j)*dt)*u(i,j,k);
                coe2 = (1-A(i,j)*dt)/(1+A(i,j)*dt)*u(i,j,k-1);
                coe3 = (a(1)*(u(i+1,j,k)-2*u(i,j,k)+u(i-1,j,k))+a(2)*(u(i+2,j,k)-2*u(i,j,k)+u(i-2,j,k))+a(3)*(u(i+3,j,k)-2*u(i,j,k)+u(i-3,j,k))+...
                    a(4)*(u(i+4,j,k)-2*u(i,j,k)+u(i-4,j,k))+a(5)*(u(i+5,j,k)-2*u(i,j,k)+u(i-5,j,k)))/(dx*dx);
                coe4 = (a(1)*(u(i,j+1,k)-2*u(i,j,k)+u(i,j-1,k))+a(2)*(u(i,j+2,k)-2*u(i,j,k)+u(i,j-2,k))+a(3)*(u(i,j+3,k)-2*u(i,j,k)+u(i,j-3,k))+...
                    a(4)*(u(i,j+4,k)-2*u(i,j,k)+u(i,j-4,k))+a(5)*(u(i,j+5,k)-2*u(i,j,k)+u(i,j-5,k)))/(dz*dz);
                u(i,j,k+1) = coe1 - coe2 + (v(i,j)*v(i,j)*dt*dt)/(1+A(i,j)*dt)*(coe3+coe4);                             
            end
        end
    end
    if j == start_traces+(traces-1)*traces_interval
        seis_wave = u(n+PML_n+1,j,k);
    end
    k
end

% u1 = u./(max(max(max(abs(u)))));

filename='��ά��������PML_2_10b_������˥������500.gif';  % ͼƬ�ļ���
for k=20:20:Nt
%     vv=imread('two_layers_velocity.png');imshow(vv);
%     hold on
    pcolor(u(n+1:sum_Nx-n,n+1:sum_Nz-n,k))
    %     pcolor(u1(n+PML_n+1:n+1+Nx+PML_n,n+PML_n+1:n+1+Nz+PML_n,k))
    %     pcolor(P_wave(:,:,k));
    shading interp;
    %     colormap('gray');
%     axis equal;
%     axis([n,sum_Nx-n,n,sum_Nz-n]);
    caxis([-1 1].*0.03)
    set(gca,'Ydir','reverse','xAxislocation','top');
    xlabel('x','fontname','times new roman','fontsize',30); ylabel('z','fontname','times new roman','fontsize',30);
    set(gca,'fontname','times new roman','fontsize',30)
    set(gcf,'unit','centimeters','position',[5 2 30 20])
        title(sprintf('time=%4f s',k*dt));
    % if(k==201)
    %     keyboard;
    % end
    f=getframe(gcf);     % ��ȡ�����������ݵ�ͼ��
    imind=frame2im(f);   %
    [imind,cm] = rgb2ind(imind,256);   %[X,MAP] =rgb2ind(RGB,N)�����RGBͼ��ת��Ϊ����ͼ��X������Nָ��MAP����ɫ����������ֵ��MAP��Ϊ����ͼ��ĵ�ɫ��
    if k==20
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.05); %��ͼ������imind���������ɫͼcmд���µ��ļ��У������ӳ�ʱ��Ϊ0.05��д��������ļ�
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.001); %�����ӳ�ʱ��Ϊ0.1��д��������ļ�
    end
end

%%
for i=1:Nt
    for j=1:traces
        seis_wave(j,i)=u(n+PML_n+1,start_traces+(j-1)*traces_interval,i);
    end
end
figure(4);
wigb(seis_wave')
xlabel('Trace','fontname','times new roman','fontsize',30); ylabel('T/ms','fontname','times new roman','fontsize',30);
set(gca,'fontname','times new roman','fontsize',30)
set(gcf,'unit','centimeters','position',[2 2 25 20])
saveas(gcf,'seis_wave.png')
%%
geo_location_x = [100:5:195];
geo_location_y = ones(1,traces);
figure(3)
pcolor(u1(n+PML_n+1:n+1+Nx+PML_n,n+PML_n+1:n+1+Nz+PML_n,601))
shading interp;
colormap('gray');
axis equal;
axis([1,Nx,1,Nz]);
set(gca,'Ydir','reverse','xAxislocation','top');
xlabel('x','fontname','times new roman','fontsize',30); ylabel('z','fontname','times new roman','fontsize',30);
set(gca,'fontname','times new roman','fontsize',30)
set(gcf,'unit','centimeters','position',[2 2 25 20])
hold on;
scatter(geo_location_x,geo_location_y,200,'rv','filled')

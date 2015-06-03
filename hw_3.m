clc
clear all
close all

format short
%%Data
mu=398600.2;
Re=6378.16;
T=100*60;
a=((T/(2*pi))^2*mu)^(1/3);
alt=a-Re;%km
n=1/sqrt(a^3/mu);
m_start=100;
Isp_chem=200;
g0=9.81;

tf=T/60;

G=[0,0,0;0,0,0;0,0,0;1,0,0;0,1,0;0,0,1];
Gtrans=G';
F=[0,0,0,1,0,0;0,0,0,0,1,0;0,0,0,0,0,1; 3*n^2, 0, 0,0,2*n,0;0,0,0,-2*n,0,0;0,0,-n*n,0,0,0];
% G2=[0,0;0,0;1,0;0,1];
%  F=[0,0,1,0;0,0,0,1;3*n^2, 0,0,2*n;0,0,-2*n,0];
syms tao t
M=-int(expm(F*(t-tao))*G*G'*expm(-F'*(tao)),tao,0,t);
Mg=matlabFunction(M);
Psi=eye(6);
%%
X0=[0.5,0,0.3,0,0,0];%initial conditions
yf=[-60/sqrt(3),-60/sqrt(3),-60/sqrt(3),0.04/sqrt(3),0.04/sqrt(3),0.04/sqrt(3)]*10^-3;%constrains
ind=0;
fsize_0=0;
dt=0.1; %time step
for times=0:dt:tf-dt  %numerical integration
Mf=real(double(Mg(tf-times)));
ind=ind+1;
K=double(G'*expm(F'*(tf-times))*...
    (inv(Mf*expm(F'*(tf-times)))));
f(:,ind)=K*(expm(F*(tf-times)))*X0'-K*yf';
Xt(:,ind)=(eye(6)+(F+G*K*expm(F*(tf-times)))*dt)*X0'-G*K*yf'*dt;
X0=Xt(:,ind)';
fsize(ind)=norm(f(:,ind));
deltaV=fsize(ind)*dt+fsize_0;
fsize_0=deltaV;
end
disp(deltaV);
figure; plot3(Xt(1,:),Xt(2,:),Xt(3,:));grid on
xlabel('X [Km]');ylabel('Y [Km]');zlabel('Z [Km]');title(['Closed loop optimization - back to the space station trajectory, t_f=',num2str(tf),' sec']);
hold on;
g=0;N=40;
for i=0:360/N:360 
  
    for j=0:360/N:360
        g=g+1;
        ball(g,:)=[ 0.06*cosd(j)*cosd(i),0.06*cosd(i)*sind(j),0.06*sind(i)]/sqrt(3);
    end;
end
  hold on
  plot3(ball(:,1),ball(:,2),ball(:,3),'b','linewidth',2);
axis equal
grid on;


figure;
subplot(3,1,1); 
plot(0:dt:(tf-dt),f(1,1:end)); grid on;
ylabel('f_x [Km/s^2]');xlabel('Time [s]');title('Closed loop optimization - back to the space station acceleration');


subplot(3,1,2); 
plot(0:dt:(tf-dt),f(2,1:end)); grid on;
ylabel('f_y [Km/s^2]');xlabel('Time [s]');

subplot(3,1,3);
plot(0:dt:(tf-dt),f(3,1:end)); grid on;
ylabel('f_z [Km/s^2]');xlabel('Time [s]');
pos1=Xt(:,end)
delv1=deltaV
errorpos=abs(pos1+[0.06/sqrt(3)*ones(3,1);-0.00004/sqrt(3)*ones(3,1)])*1000;
index=ind;
Xerror=norm(errorpos(1:3))
Verror=norm(errorpos(4:6))
%% from safety zone edge to [0,0,0]

t02=tf-dt; %time gap
tf=tf+T/4; %Time for strait line with onstant speed R/V=t
X0=Xt(:,end)';%initial conditions
yf=zeros(1,6); %target point
for times=t02:dt:tf-dt
Mf=real(double(Mg(tf-times)));
ind=ind+1;
K=double(G'*expm(F'*(tf-times))*...
    (inv(Mf*expm(F'*(tf-times)))));
f(:,ind)=K*(expm(F*(tf-times)))*X0'-K*yf';
Xt(:,ind)=(eye(6)+(F+G*K*(expm(F*(tf-times))))*dt)*X0'-G*K*yf'*dt;
X0=Xt(:,ind)';
fsize(ind)=norm(f(:,ind));
deltaV=fsize(ind)*dt+fsize_0;
fsize_0=deltaV;
end

figure; plot3(Xt(1,:),Xt(2,:),Xt(3,:));grid on
 xlabel('X [Km]');ylabel('Y [Km]');zlabel('Z [Km]');title(['Closed loop optimization - back to the space station trajectory, ',num2str(tf),' sec']);

 figure; plot3(Xt(1,index:end),Xt(2,index:end),Xt(3,index:end));grid on
 xlabel('X [Km]');ylabel('Y [Km]');zlabel('Z [Km]');title('Closed loop optimization - back to the space station trajectory');
 % hold on;
% g=0;N=40;
%   
%     for j=0:360/N:360
%         g=g+1;
%         ball(g,:)=[ 0.06*cosd(j)*cosd(i),0.06*cosd(i)*sind(j),0.06*sind(i)]/sqrt(3);
%     end;
% end
%  hold on
%   plot3(ball(:,1),ball(:,2),ball(:,3),'b','linewidth',2);
% axis equal
% grid on;

maxf=max(fsize)
delv2=deltaV

figure; 
subplot(3,1,1); 
plot(0:dt:(tf),f(1,1:end));grid on;
ylabel('f_x [Km/s^2]');xlabel('Time [s]');title('Closed loop optimization - back to the space station acceleration');

subplot(3,1,2); 
plot(0:dt:(tf),f(2,1:end));grid on;
ylabel('f_y [Km/s^2]');xlabel('Time [s]');

subplot(3,1,3); 
plot(0:dt:(tf),f(3,1:end));grid on;
ylabel('f_z [Km/s^2]');xlabel('Time [s]');

X2error=norm(Xt(1:3,end))
V2error=norm(Xt(4:6,end))
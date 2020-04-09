clear all, close all, clc

m = 1;
M = 5;
L = 2;
g = 10;
d = 1;
I = 4/3*m*L^2;

xf = [0; 0; pi; 0];  % fixed point to linearize around
uf = 0;

% A = [0 1 0 0;
%     0 d*(I + L^2*m)/(L^2*m^2*cos(xf(3))^2 - (I + L^2*m)*(M + m)) -L*m*(L*m*(L^2*g*m^2*sin(2*xf(3)) + 2*(I + L^2*m)*(L*m*xf(4)^2*sin(xf(3)) - d*xf(2) + uf))*sin(xf(3))*cos(xf(3)) + (L^2*m^2*cos(xf(3))^2 - (I + L^2*m)*(M + m))*(L*g*m*cos(2*xf(3)) + xf(4)^2*(I + L^2*m)*cos(xf(3))))/(L^2*m^2*cos(xf(3))^2 - (I + L^2*m)*(M + m))^2 -2*L*m*xf(4)*(I + L^2*m)*sin(xf(3))/(L^2*m^2*cos(xf(3))^2 - (I + L^2*m)*(M + m));
%     0 0 0 1;
%     0 L*d*m*cos(xf(3))/(I*M + I*m + L^2*M*m + L^2*m^2*sin(xf(3))^2) L*m*(L^2*m^2*(L*m*xf(4)^2*sin(2*xf(3)) + 2*M*g*sin(xf(3)) - 2*d*xf(2)*cos(xf(3)) + 2*g*m*sin(xf(3)) + 2*uf*cos(xf(3)))*sin(xf(3))*cos(xf(3)) - (I*M + I*m + L^2*M*m + L^2*m^2*sin(xf(3))^2)*(L*m*xf(4)^2*cos(2*xf(3)) + M*g*cos(xf(3)) + d*xf(2)*sin(xf(3)) + g*m*cos(xf(3)) - uf*sin(xf(3))))/(I*M + I*m + L^2*M*m + L^2*m^2*sin(xf(3))^2)^2 -L^2*m^2*xf(4)*sin(2*xf(3))/(I*M + I*m + L^2*M*m + L^2*m^2*sin(xf(3))^2)];
% B = [0; -(I + L^2*m)/(L^2*m^2*cos(xf(3))^2 - (I + L^2*m)*(M + m)); 0; -L*m*cos(xf(3))/(I*M + I*m + L^2*M*m + L^2*m^2*sin(xf(3))^2)];

A = [0 1 0;
    L*m*(L^2*m^2*(L*m*xf(4)^2*sin(2*xf(3)) + 2*M*g*sin(xf(3)) - 2*d*xf(2)*cos(xf(3)) + 2*g*m*sin(xf(3)) + 2*uf*cos(xf(3)))*sin(xf(3))*cos(xf(3)) - (I*M + I*m + L^2*M*m + L^2*m^2*sin(xf(3))^2)*(L*m*xf(4)^2*cos(2*xf(3)) + M*g*cos(xf(3)) + d*xf(2)*sin(xf(3)) + g*m*cos(xf(3)) - uf*sin(xf(3))))/(I*M + I*m + L^2*M*m + L^2*m^2*sin(xf(3))^2)^2 -L^2*m^2*xf(4)*sin(2*xf(3))/(I*M + I*m + L^2*M*m + L^2*m^2*sin(xf(3))^2) 0;
    1 0 0];
B = [0; -L*m*cos(xf(3))/(I*M + I*m + L^2*M*m + L^2*m^2*sin(xf(3))^2); 0];
Afull = [0 1 0 0 0;
    0 d*(I + L^2*m)/(L^2*m^2*cos(xf(3))^2 - (I + L^2*m)*(M + m)) -L*m*(L*m*(L^2*g*m^2*sin(2*xf(3)) + 2*(I + L^2*m)*(L*m*xf(4)^2*sin(xf(3)) - d*xf(2) + uf))*sin(xf(3))*cos(xf(3)) + (L^2*m^2*cos(xf(3))^2 - (I + L^2*m)*(M + m))*(L*g*m*cos(2*xf(3)) + xf(4)^2*(I + L^2*m)*cos(xf(3))))/(L^2*m^2*cos(xf(3))^2 - (I + L^2*m)*(M + m))^2 -2*L*m*xf(4)*(I + L^2*m)*sin(xf(3))/(L^2*m^2*cos(xf(3))^2 - (I + L^2*m)*(M + m)) 0;
    0 0 0 1 0;
    0 L*d*m*cos(xf(3))/(I*M + I*m + L^2*M*m + L^2*m^2*sin(xf(3))^2) L*m*(L^2*m^2*(L*m*xf(4)^2*sin(2*xf(3)) + 2*M*g*sin(xf(3)) - 2*d*xf(2)*cos(xf(3)) + 2*g*m*sin(xf(3)) + 2*uf*cos(xf(3)))*sin(xf(3))*cos(xf(3)) - (I*M + I*m + L^2*M*m + L^2*m^2*sin(xf(3))^2)*(L*m*xf(4)^2*cos(2*xf(3)) + M*g*cos(xf(3)) + d*xf(2)*sin(xf(3)) + g*m*cos(xf(3)) - uf*sin(xf(3))))/(I*M + I*m + L^2*M*m + L^2*m^2*sin(xf(3))^2)^2 -L^2*m^2*xf(4)*sin(2*xf(3))/(I*M + I*m + L^2*M*m + L^2*m^2*sin(xf(3))^2) 0;
    0 0 1 0 0];
Bfull = [0; -(I + L^2*m)/(L^2*m^2*cos(xf(3))^2 - (I + L^2*m)*(M + m)); 0; -L*m*cos(xf(3))/(I*M + I*m + L^2*M*m + L^2*m^2*sin(xf(3))^2); 0];


eig(A)
rank(ctrb(A,B))

%%  Design LQR controller
% Q = [1e-16 0 0 0;
%     0 1 0 0;
%     0 0 1 0;
%     0 0 0 1];
% R = 1e-4;

Q = [1 0 0;
    0 1 0;
    0 0 100];
R = 1e-3;

K = lqr(A,B,Q,R);

%% Simulate closed-loop system
tspan = 0:.001:10;
% x0 = [-1; 0; pi+0.7; 0];  % initial condition
% wr = [1; 0; pi; 0];      % reference position

x0 = [-1; 0; pi+0.7; 0; 0];  % initial condition
wr = [1; 0; pi; 0; 0];      % reference position
% u=@(x)-K*(reshape(x,length(wr),1) - wr);       % control law

u=@(x)-K*(reshape(x(3:end),length(wr(3:end)),1) - wr(3:end)); % control law
% [t,x] = ode45(@(t,x)pendcart(x,m,M,L,g,d,u(x)),tspan,x0);

[t,x] = ode45(@(t,x)pendcart(x,m,M,L,g,d,u(x),wr),tspan,x0);

for k=1:100:length(t)
    drawpend(x(k,:),m,M,L);
end

%%
plot(t,x,'LineWidth',2); hold on
% l1 = legend('x','v','\theta','\omega');
l1 = legend('x','v','\theta','\omega','i');
set(l1,'Location','SouthEast')
set(gcf,'Position',[100 100 500 200])
xlabel('Time')
ylabel('State')
grid on
set(gcf,'PaperPositionMode','auto')
% print('-depsc2', '-loose', 'figures/FIG_02_LQR');

%% Compare with many examples of Pole Placement
% K = lqr(A,B,Q,R);
% u=@(x)-K*(x - wr);       % control law
% [t,x] = ode45(@(t,x)pendcart(x,m,M,L,g,d,u(x)),tspan,x0);
% xLQR = x;
% for k=1:length(t)
%     JLQR(k) = (x(k,:)-wr')*Q*(x(k,:)'-wr) + u(x(k,:)')^2*R;
% end
% 
% CC = [0    0.4470    0.7410
%     0.8500    0.3250    0.0980
%     0.9290    0.6940    0.1250
%     0.4940    0.1840    0.5560
%     0.4660    0.6740    0.1880
%     0.3010    0.7450    0.9330
%     0.6350    0.0780    0.1840];
% 
% CCgray = [0.2    0.6470    0.9410
%     0.9500    0.4250    0.1980
%     1    0.7940    0.2250
%     0.5940    0.2840    0.6560
%     0.4660    0.6740    0.1880
%     0.3010    0.7450    0.9330
%     0.6350    0.0780    0.1840];
% 
% for count = 1:100
%     p = [-.5-3*rand; -.5-3*rand; -.5-3*rand; -.5-3*rand];
% 	K = place(A,B,p);
%     u=@(x)-K*(x - wr);       % control law
%     [t,x] = ode45(@(t,x)pendcart(x,m,M,L,g,d,u(x)),tspan,x0);
%     figure(1)
%     for j=1:4
%         plot(t(1:50:end),x(1:50:end,j),'Color',[.2 .2 .2] + .5*CC(j,:)), hold on;
%     end
%     for k=1:length(t)
%         J(k) = (x(k,:)-wr')*Q*(x(k,:)'-wr) + u(x(k,:)')^2*R;
%     end
%     figure(2)
%     Jz = cumtrapz(t,J);
%     plot(t(1:50:end),Jz(1:50:end),'Color',[.5 .5 .5]), hold on;
% end
% figure(1)
%     for j=1:4
%         plot(t(1:10:end),xLQR(1:10:end,j),'Color',CC(j,:),'LineWidth',2)
%     end
% figure(2)
% plot(t,cumtrapz(t,JLQR),'k','LineWidth',2)
% 
% figure(1)
% set(gcf,'Position',[100 100 500 300])
% xlabel('Time')
% ylabel('State')
% grid on
% set(gcf,'PaperPositionMode','auto')
% print('-depsc2', '-loose', 'figures/FIG_02b_LQRb');
% 
% figure(2)
% set(gcf,'Position',[100 100 500 300])
% xlabel('Time')
% ylabel('Cost')
% grid on
% set(gcf,'PaperPositionMode','auto')
% print('-depsc2', '-loose', 'figures/FIG_02c_LQRb');
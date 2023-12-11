%% 
close all
clear
clc

% ut+x*ux=t
% u(x0,0) = u0(x) = sin(2*x)

Nx = 21;
Nt = 51;
x0 = linspace(-1,1,Nx)';
t  = linspace(0,1,Nt)';
for j=1:length(t)
    x(:,j) = x0*exp(t(j));
    u(:,j) = sin(2*x0) + t(j)^2/2;
end

subplot(131)
plot(x(:,1),u(:,1)); hold on
plot(x(:,31),u(:,31)); hold on
plot(x(:,51),u(:,51)); hold on
xlabel('x')
ylabel('u(x,t)')
legend('t_0','t_1','t_2','location','best')
subplot(132)
for j=1:length(t)
    plot3(t(j)*ones(Nx,1),x(:,j),u(:,j)); hold on
end
for i=1:length(x0)
    plot3(t,x(i,:),-1*ones(Nt,1),'k'); hold on
end
xlabel('t')
ylabel('x')
zlabel('u(x,t)')
subplot(133)
for i=1:length(x0)
    plot3(t,x(i,:),u(i,:),'k'); hold on
    plot3(t,x(i,:),-1*ones(Nt,1),'k'); hold on
end
xlabel('t')
ylabel('x')
zlabel('u(x,t)')

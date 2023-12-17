%% ========================================================================
close all
clear
clc
setup
label    = 'sine';
iFig     = 1;
saveFigs = 1;

% PARAMETERS --------------------------------------------------------------
timeMethod   = 3;
flowType     = 4;
dataFlow.u0  = 1;
dataFlow.k   = 2;
dataFlow.A   = 0.5;
dataFlow.w   = 0;
dataFlow.phi = 0;
f1Type     = 3; % Stokes random
dataf1     = [];
nt         = [301 10];
tLim       = [0 10];

% Forcing and initial condition
mean_a1  = 1;
sigma_a1 = 0.1;
a1a1     = sigma_a1^2;
cm3_a1   = 0;
lim1     = mean_a1-0.5*sqrt(12)*sigma_a1;
lim2     = mean_a1+0.5*sqrt(12)*sigma_a1;
cm4_a1   = (lim2-lim1)^4/80;
xp00     = 0;
yp00     = 0;
up00     = 0;
vp00     = 0;


%% MC-PSIC ================================================================

% -------------------------------------------------------------------------
taup = 1;
ns   = 10000; % Numer of samples per cloud in MC-PSIC
a1 = linspace(mean_a1-0.5*sqrt(12)*sigma_a1,mean_a1+0.5*sqrt(12)*sigma_a1,ns)';
for k=1:ns
    xp0(1,1,k)       = xp00;
    yp0(1,1,k)       = yp00;
    up0(1,1,k)       = up00;
    vp0(1,1,k)       = vp00;
    dataf1.a1(1,1,k) = a1(k);
end
M1 = FUNC_solverMCPSIC2D_inertial(xp0,yp0,up0,vp0,tLim,nt,flowType,dataFlow,f1Type,dataf1,timeMethod,taup);


%% MoP ====================================================================

taup = 1;
clear xp0 yp0 up0 vp0
xp0(1,1) = xp00;
yp0(1,1) = yp00;
up0(1,1) = up00;
vp0(1,1) = vp00;

% M=1 ---------------------------------------------------------------------
dataf1.mean_a1  = mean_a1;
dataf1.a1       = a1;
dataf1.a1a1     = a1a1;
dataf1.a1a1a1   = cm3_a1;
dataf1.a1a1a1a1 = cm4_a1*0;
P1 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,tLim,nt, ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

% M=2 ---------------------------------------------------------------------
M = 2; k = 1;
for ia=1:M
    if M == 1
        La = 0;
        Da = 0;
    elseif M == 2
        La = 2*sqrt( 0.5*a1a1 );
        Da = 0.5*La;
    elseif M == 3
        La = sqrt( a1a1 );
        Da = La;
    elseif M == 4
        La = sqrt( 3*a1a1/(2*1.5^2+2*0.5^2) );
        Da = 1.5*La;
    elseif M == 5
        La = sqrt( 4/10 *a1a1 );
        Da = 2*La;
    elseif M == 6
        La = sqrt( 5/2 *a1a1/(2.5^2+1.5^2+0.5^2) );
        Da = 2.5*La;
    elseif M == 7
        La = sqrt( 6/2 *a1a1/(3^2+2^2+1^2) );
        Da = 3*La;
    elseif M == 8
        La = sqrt( 7/2 *a1a1/(3.5^2+2.5^2+1.5^2+0.5^2) );
        Da = 3.5*La;
    end
    mean_a_s(k,1) = mean_a1 - Da + (ia-1)*La;
    a1a1_s(k,1)   = a1a1/M;
    k = k + 1;
end
dataf1.mean_a1  = mean_a_s(1);
dataf1.a1       = a1;
dataf1.a1a1     = a1a1_s(1);
dataf1.a1a1a1   = cm3_a1;
dataf1.a1a1a1a1 = cm4_a1*0;
P21 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,tLim,nt, ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

dataf1.mean_a1  = mean_a_s(2);
dataf1.a1       = a1;
dataf1.a1a1     = a1a1_s(2);
dataf1.a1a1a1   = cm3_a1;
dataf1.a1a1a1a1 = cm4_a1*0;
P22 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,tLim,nt, ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

% Joining 
w = 1/M; % the weights would be all equal
for j=1:nt(1)
    PP2.mean_xp(j,1) = w*P21.MoM.mean_xp(1,1,j)+w*P22.MoM.mean_xp(1,1,j);
    PP2.mean_up(j,1) = w*P21.MoM.mean_up(1,1,j)+w*P22.MoM.mean_up(1,1,j);
    PP2.xpxp(j,1)    = w*P21.MoM.xpxp(1,1,j)+w*P22.MoM.xpxp(1,1,j) ...
        +w*(P21.MoM.mean_xp(1,1,j)-PP2.mean_xp(j)).^2+w*(P22.MoM.mean_xp(1,1,j)-PP2.mean_xp(j)).^2;
    PP2.upup(j,1)    = w*P21.MoM.upup(1,1,j)+w*P22.MoM.upup(1,1,j) ...
        +w*(P21.MoM.mean_up(1,1,j)-PP2.mean_up(j)).^2+w*(P22.MoM.mean_up(1,1,j)-PP2.mean_up(j)).^2;
end

% M=3 ---------------------------------------------------------------------
M = 3; k = 1;
for ia=1:M
    if M == 1
        La = 0;
        Da = 0;
    elseif M == 2
        La = 2*sqrt( 0.5*a1a1 );
        Da = 0.5*La;
    elseif M == 3
        La = sqrt( a1a1 );
        Da = La;
    elseif M == 4
        La = sqrt( 3*a1a1/(2*1.5^2+2*0.5^2) );
        Da = 1.5*La;
    elseif M == 5
        La = sqrt( 4/10 *a1a1 );
        Da = 2*La;
    elseif M == 6
        La = sqrt( 5/2 *a1a1/(2.5^2+1.5^2+0.5^2) );
        Da = 2.5*La;
    elseif M == 7
        La = sqrt( 6/2 *a1a1/(3^2+2^2+1^2) );
        Da = 3*La;
    elseif M == 8
        La = sqrt( 7/2 *a1a1/(3.5^2+2.5^2+1.5^2+0.5^2) );
        Da = 3.5*La;
    end
    mean_a_s(k,1) = mean_a1 - Da + (ia-1)*La;
    a1a1_s(k,1)   = a1a1/M;
    k = k + 1;
end
dataf1.mean_a1  = mean_a_s(1);
dataf1.a1       = a1;
dataf1.a1a1     = a1a1_s(1);
dataf1.a1a1a1   = cm3_a1;
dataf1.a1a1a1a1 = cm4_a1*0;
P31 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,tLim,nt, ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

dataf1.mean_a1  = mean_a_s(2);
dataf1.a1       = a1;
dataf1.a1a1     = a1a1_s(2);
dataf1.a1a1a1   = cm3_a1;
dataf1.a1a1a1a1 = cm4_a1*0;
P32 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,tLim,nt, ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

dataf1.mean_a1  = mean_a_s(3);
dataf1.a1       = a1;
dataf1.a1a1     = a1a1_s(3);
dataf1.a1a1a1   = cm3_a1;
dataf1.a1a1a1a1 = cm4_a1*0;
P33 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,tLim,nt, ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

% Joining 
w = 1/M; % the weights would be all equal
for j=1:nt(1)
    PP3.mean_xp(j,1) = w*P31.MoM.mean_xp(1,1,j)+w*P32.MoM.mean_xp(1,1,j)+w*P33.MoM.mean_xp(1,1,j);
    PP3.mean_up(j,1) = w*P31.MoM.mean_up(1,1,j)+w*P32.MoM.mean_up(1,1,j)+w*P33.MoM.mean_up(1,1,j);
    PP3.xpxp(j,1)    = w*P31.MoM.xpxp(1,1,j)+w*P32.MoM.xpxp(1,1,j)+w*P33.MoM.xpxp(1,1,j) ...
        +w*(P31.MoM.mean_xp(1,1,j)-PP3.mean_xp(j)).^2+w*(P32.MoM.mean_xp(1,1,j)-PP3.mean_xp(j)).^2+w*(P33.MoM.mean_xp(1,1,j)-PP3.mean_xp(j)).^2;
    PP3.upup(j,1)    = w*P31.MoM.upup(1,1,j)+w*P32.MoM.upup(1,1,j)+w*P33.MoM.upup(1,1,j) ...
        +w*(P31.MoM.mean_up(1,1,j)-PP3.mean_up(j)).^2+w*(P32.MoM.mean_up(1,1,j)-PP3.mean_up(j)).^2+w*(P33.MoM.mean_up(1,1,j)-PP3.mean_up(j)).^2;
end

% PLOTS -------------------------------------------------------------------

close all
plot(squeeze(M1.upup)); hold on
plot(squeeze(P1.MoM.upup),'--'); hold on
plot(PP2.upup,':'); hold on

close all
subplot(221)
plot(squeeze(M1.mean_xp)); hold on
plot(squeeze(P1.MoM.mean_xp),'--'); hold on
plot(PP2.mean_xp,'-.'); hold on
plot(PP3.mean_xp,':'); hold on
subplot(222)
plot(squeeze(M1.mean_up)); hold on
plot(squeeze(P1.MoM.mean_up),'--'); hold on
plot(PP2.mean_up,'-.'); hold on
plot(PP3.mean_up,':'); hold on
subplot(223)
plot(squeeze(M1.xpxp)); hold on
plot(squeeze(P1.MoM.xpxp),'--'); hold on
plot(PP2.xpxp,'-.'); hold on
plot(PP3.xpxp,':'); hold on
subplot(224)
plot(squeeze(M1.upup)); hold on
plot(squeeze(P1.MoM.upup),'--'); hold on
plot(PP2.upup,'-.'); hold on
plot(PP3.upup,':'); hold on


e1 = norm(squeeze(P1.MoM.upup.^0.5)-squeeze(M1.upup.^0.5),2)/norm(squeeze(M1.upup.^0.5),Inf);
e2 = norm(squeeze(PP2.upup.^0.5)-squeeze(M1.upup.^0.5),2)/norm(squeeze(M1.upup.^0.5),Inf);
e3 = norm(squeeze(PP3.upup.^0.5)-squeeze(M1.upup.^0.5),2)/norm(squeeze(M1.upup.^0.5),Inf);

E1 = norm(squeeze(P1.MoM.upup.^0.5)-PP3.upup.^0.5,2)/norm(PP3.upup.^0.5,Inf);
E2 = norm(squeeze(PP2.upup.^0.5)   -PP3.upup.^0.5,2)/norm(PP3.upup.^0.5,Inf);


x = linspace(1,3,100);
y = 0.8*x.^(-2);
close all
subplot(121)
loglog([1 2 3],[e1 e2 e3],'o--'); hold on
loglog(x,y,'--k')
subplot(122)
loglog([1 2],[E1 E2],'o--'); hold on
loglog(x,y,'--k')




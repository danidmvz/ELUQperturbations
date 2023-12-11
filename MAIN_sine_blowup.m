%% ========================================================================
close all
clear
clc
setup
label    = 'sine';
iFig     = 1;
saveFigs = 0;

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
nt         = [301 2];
tLim       = [0 20];
taup       = 1;

% Forcing and initial condition
sigma_a1 = 0.1;
cm3_a1   = 0;
cm4_a1   = 0;
mean_a1  = 1;
xp00     = 0;
yp00     = 0;
up00     = 0;
vp00     = 0;


%% MC-PSIC ================================================================

% -------------------------------------------------------------------------
taup = 1;
ns   = 10000; % Numer of samples per cloud in MC-PSIC
a1   = linspace(mean_a1-0.5*sqrt(12)*sigma_a1,mean_a1+0.5*sqrt(12)*sigma_a1,ns)';
for k=1:ns
    xp0(1,1,k)       = xp00;
    yp0(1,1,k)       = yp00;
    up0(1,1,k)       = up00;
    vp0(1,1,k)       = vp00;
    dataf1.a1(1,1,k) = a1(k);
end
M1 = FUNC_solverMCPSIC2D_inertial(xp0,yp0,up0,vp0,[tLim(1) taup*tLim(2)],[nt 2],flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

% -------------------------------------------------------------------------
taup = 2;
ns   = 10000; % Numer of samples per cloud in MC-PSIC
a1   = linspace(mean_a1-0.5*sqrt(12)*sigma_a1,mean_a1+0.5*sqrt(12)*sigma_a1,ns)';
for k=1:ns
    xp0(1,1,k)       = xp00;
    yp0(1,1,k)       = yp00;
    up0(1,1,k)       = up00;
    vp0(1,1,k)       = vp00;
    dataf1.a1(1,1,k) = a1(k);
end
M2 = FUNC_solverMCPSIC2D_inertial(xp0,yp0,up0,vp0,[tLim(1) taup*tLim(2)],[taup*nt 2],flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

% -------------------------------------------------------------------------
taup = 3;
ns   = 10000; % Numer of samples per cloud in MC-PSIC
a1   = linspace(mean_a1-0.5*sqrt(12)*sigma_a1,mean_a1+0.5*sqrt(12)*sigma_a1,ns)';
for k=1:ns
    xp0(1,1,k)       = xp00;
    yp0(1,1,k)       = yp00;
    up0(1,1,k)       = up00;
    vp0(1,1,k)       = vp00;
    dataf1.a1(1,1,k) = a1(k);
end
M3 = FUNC_solverMCPSIC2D_inertial(xp0,yp0,up0,vp0,[tLim(1) taup*tLim(2)],[taup*nt 2],flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

% -------------------------------------------------------------------------
taup = 4;
ns   = 10000; % Numer of samples per cloud in MC-PSIC
a1   = linspace(mean_a1-0.5*sqrt(12)*sigma_a1,mean_a1+0.5*sqrt(12)*sigma_a1,ns)';
for k=1:ns
    xp0(1,1,k)       = xp00;
    yp0(1,1,k)       = yp00;
    up0(1,1,k)       = up00;
    vp0(1,1,k)       = vp00;
    dataf1.a1(1,1,k) = a1(k);
end
M4 = FUNC_solverMCPSIC2D_inertial(xp0,yp0,up0,vp0,[tLim(1) taup*tLim(2)],[taup*nt 2],flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

% -------------------------------------------------------------------------
taup = 5;
ns   = 10000; % Numer of samples per cloud in MC-PSIC
a1   = linspace(mean_a1-0.5*sqrt(12)*sigma_a1,mean_a1+0.5*sqrt(12)*sigma_a1,ns)';
for k=1:ns
    xp0(1,1,k)       = xp00;
    yp0(1,1,k)       = yp00;
    up0(1,1,k)       = up00;
    vp0(1,1,k)       = vp00;
    dataf1.a1(1,1,k) = a1(k);
end
M5 = FUNC_solverMCPSIC2D_inertial(xp0,yp0,up0,vp0,[tLim(1) taup*tLim(2)],[taup*nt 2],flowType,dataFlow,f1Type,dataf1,timeMethod,taup);



%% MoP ====================================================================

% -------------------------------------------------------------------------
taup = 1;
clear xp0 yp0 up0 vp0
xp0(1,1) = xp00;
yp0(1,1) = yp00;
up0(1,1) = up00;
vp0(1,1) = vp00;
dataf1.a1       = a1;
dataf1.a1a1     = sigma_a1^2;
dataf1.a1a1a1   = cm3_a1;
dataf1.a1a1a1a1 = cm4_a1;
P1 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,[tLim(1) taup*tLim(2)],[taup*nt 2], ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

% -------------------------------------------------------------------------
taup = 2;
clear xp0 yp0 up0 vp0
xp0(1,1) = xp00;
yp0(1,1) = yp00;
up0(1,1) = up00;
vp0(1,1) = vp00;
dataf1.a1       = a1;
dataf1.a1a1     = sigma_a1^2;
dataf1.a1a1a1   = cm3_a1;
dataf1.a1a1a1a1 = cm4_a1;
P2 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,[tLim(1) taup*tLim(2)],[taup*nt 2], ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

% -------------------------------------------------------------------------
% taup = 2.5;
% clear xp0 yp0 up0 vp0
% xp0(1,1) = xp00;
% yp0(1,1) = yp00;
% up0(1,1) = up00;
% vp0(1,1) = vp00;
% dataf1.a1       = a1;
% dataf1.a1a1     = sigma_a1^2;
% dataf1.a1a1a1   = cm3_a1;
% dataf1.a1a1a1a1 = cm4_a1;
% P25 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,tLim,nt, ...
%     flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

% -------------------------------------------------------------------------
taup = 3;
clear xp0 yp0 up0 vp0
xp0(1,1) = xp00;
yp0(1,1) = yp00;
up0(1,1) = up00;
vp0(1,1) = vp00;
dataf1.a1       = a1;
dataf1.a1a1     = sigma_a1^2;
dataf1.a1a1a1   = cm3_a1;
dataf1.a1a1a1a1 = cm4_a1;
P3 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,[tLim(1) taup*tLim(2)],[taup*nt 2], ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

% -------------------------------------------------------------------------
% taup = 3.5;
% clear xp0 yp0 up0 vp0
% xp0(1,1) = xp00;
% yp0(1,1) = yp00;
% up0(1,1) = up00;
% vp0(1,1) = vp00;
% dataf1.a1       = a1;
% dataf1.a1a1     = sigma_a1^2;
% dataf1.a1a1a1   = cm3_a1;
% dataf1.a1a1a1a1 = cm4_a1;
% P35 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,tLim,nt, ...
%     flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

% -------------------------------------------------------------------------
taup = 4;
clear xp0 yp0 up0 vp0
xp0(1,1) = xp00;
yp0(1,1) = yp00;
up0(1,1) = up00;
vp0(1,1) = vp00;
dataf1.a1       = a1;
dataf1.a1a1     = sigma_a1^2;
dataf1.a1a1a1   = cm3_a1;
dataf1.a1a1a1a1 = cm4_a1;
P4 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,[tLim(1) taup*tLim(2)],[taup*nt 2], ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

% -------------------------------------------------------------------------
taup = 5;
clear xp0 yp0 up0 vp0
xp0(1,1) = xp00;
yp0(1,1) = yp00;
up0(1,1) = up00;
vp0(1,1) = vp00;
dataf1.a1       = a1;
dataf1.a1a1     = sigma_a1^2;
dataf1.a1a1a1   = cm3_a1;
dataf1.a1a1a1a1 = cm4_a1;
P5 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,[tLim(1) taup*tLim(2)],[taup*nt 2], ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

%% Orders -----------------------------------------------------------------

close all
i = 1;
j = 1;

subplot(325)
plot(P1.t,squeeze(P1.o2.xp(i,j,:)),'-','linewidth',lw,'color',color1); hold on
plot(P2.t,squeeze(P2.o2.xp(i,j,:)),'-','linewidth',lw,'color',color2); hold on
plot(P3.t,squeeze(P3.o2.xp(i,j,:)),'-','linewidth',lw,'color',color3); hold on
plot(P4.t,squeeze(P4.o2.xp(i,j,:)),'-','linewidth',lw,'color',color4); hold on
plot(P5.t,squeeze(P5.o2.xp(i,j,:)),'-','linewidth',lw,'color',color5); hold on
xlabel('t')
ylabel('O(\epsilon^2), x_p')
subplot(326)
plot(P1.t,squeeze(P1.o2.up(i,j,:)),'-','linewidth',lw,'color',color1); hold on
plot(P2.t,squeeze(P2.o2.up(i,j,:)),'-','linewidth',lw,'color',color2); hold on
plot(P3.t,squeeze(P3.o2.up(i,j,:)),'-','linewidth',lw,'color',color3); hold on
plot(P4.t,squeeze(P4.o2.up(i,j,:)),'-','linewidth',lw,'color',color4); hold on
plot(P5.t,squeeze(P5.o2.up(i,j,:)),'-','linewidth',lw,'color',color5); hold on
xlabel('t')
ylabel('O(\epsilon^2), u_p')

subplot(323)
plot(P1.t,squeeze(P1.o1.xp(i,j,:)),'-','linewidth',lw,'color',color1); hold on
plot(P2.t,squeeze(P2.o1.xp(i,j,:)),'-','linewidth',lw,'color',color2); hold on
plot(P3.t,squeeze(P3.o1.xp(i,j,:)),'-','linewidth',lw,'color',color4); hold on
plot(P4.t,squeeze(P4.o2.xp(i,j,:)),'-','linewidth',lw,'color',color5); hold on
plot(P5.t,squeeze(P5.o2.xp(i,j,:)),'-','linewidth',lw,'color',color6); hold on
xlabel('t')
ylabel('O(\epsilon), x_p')
subplot(324)
plot(P1.t,squeeze(P1.o1.up(i,j,:)),'-','linewidth',lw,'color',color1); hold on
plot(P2.t,squeeze(P2.o1.up(i,j,:)),'-','linewidth',lw,'color',color2); hold on
plot(P3.t,squeeze(P3.o1.up(i,j,:)),'-','linewidth',lw,'color',color4); hold on
plot(P4.t,squeeze(P4.o2.up(i,j,:)),'-','linewidth',lw,'color',color5); hold on
plot(P5.t,squeeze(P5.o2.up(i,j,:)),'-','linewidth',lw,'color',color6); hold on
xlabel('t')
ylabel('O(\epsilon), u_p')

subplot(321)
plot(P1.t,squeeze(P1.o0.xp(i,j,:)) ,'-','linewidth',lw,'color',color1); hold on
plot(P2.t,squeeze(P2.o0.xp(i,j,:)) ,'-','linewidth',lw,'color',color2); hold on
plot(P3.t,squeeze(P3.o0.xp(i,j,:)) ,'-','linewidth',lw,'color',color4); hold on
plot(P4.t,squeeze(P4.o2.xp(i,j,:)), ':','linewidth',lw,'color',color5); hold on
plot(P5.t,squeeze(P5.o2.xp(i,j,:)), ':','linewidth',lw,'color',color6); hold on
xlabel('t')
ylabel('O(1), x_p')
subplot(322)
plot(P1.t,squeeze(P1.o0.up(i,j,:)) ,'-','linewidth',lw,'color',color1); hold on
plot(P2.t,squeeze(P2.o0.up(i,j,:)) ,'-','linewidth',lw,'color',color2); hold on
plot(P3.t,squeeze(P3.o0.up(i,j,:)) ,'-','linewidth',lw,'color',color4); hold on
plot(P4.t,squeeze(P4.o2.up(i,j,:)), ':','linewidth',lw,'color',color5); hold on
plot(P5.t,squeeze(P5.o2.up(i,j,:)), ':','linewidth',lw,'color',color6); hold on
xlabel('t')
ylabel('O(1), u_p')



%%

close all

subplot(121)
plot(P1.t,squeeze(M1.xpxp(i,j,:).^0.5)    ,'--','linewidth',lw,'color',color14); hold on
plot(P1.t,squeeze(P1.MoM.xpxp(i,j,:).^0.5),'-.','linewidth',lw,'color',color1); hold on
plot(P2.t/2,squeeze(M2.xpxp(i,j,:).^0.5)    ,'--','linewidth',lw,'color',color14); hold on
plot(P2.t/2,squeeze(P2.MoM.xpxp(i,j,:).^0.5),'-.','linewidth',lw,'color',color2); hold on
plot(P3.t/3,squeeze(M3.xpxp(i,j,:).^0.5)    ,'--','linewidth',lw,'color',color14); hold on
plot(P3.t/3,squeeze(P3.MoM.xpxp(i,j,:).^0.5),'-.','linewidth',lw,'color',color3); hold on
plot(P4.t/4,squeeze(M4.xpxp(i,j,:).^0.5)    ,'--','linewidth',lw,'color',color14); hold on
plot(P4.t/4,squeeze(P4.MoM.xpxp(i,j,:).^0.5),'-.','linewidth',lw,'color',color4); hold on
plot(P5.t/5,squeeze(M5.xpxp(i,j,:).^0.5)    ,'--','linewidth',lw,'color',color14); hold on
plot(P5.t/5,squeeze(P5.MoM.xpxp(i,j,:).^0.5),'-.','linewidth',lw,'color',color5); hold on
subplot(122)
% close all
plot(P1.t,squeeze(M1.upup(i,j,:).^0.5)    ,'--','linewidth',lw,'color',color14); hold on
plot(P1.t,squeeze(P1.MoM.upup(i,j,:).^0.5),'-.','linewidth',lw,'color',color1); hold on
plot(P2.t/2,squeeze(M2.upup(i,j,:).^0.5)    ,'--','linewidth',lw,'color',color14); hold on
plot(P2.t/2,squeeze(P2.MoM.upup(i,j,:).^0.5),'-.','linewidth',lw,'color',color2); hold on
plot(P3.t/3,squeeze(M3.upup(i,j,:).^0.5)    ,'--','linewidth',lw,'color',color14); hold on
plot(P3.t/3,squeeze(P3.MoM.upup(i,j,:).^0.5),'-.','linewidth',lw,'color',color3); hold on
plot(P4.t/4,squeeze(M4.upup(i,j,:).^0.5)    ,'--','linewidth',lw,'color',color14); hold on
plot(P4.t/4,squeeze(P4.MoM.upup(i,j,:).^0.5),'-.','linewidth',lw,'color',color4); hold on
plot(P5.t/5,squeeze(M5.upup(i,j,:).^0.5)    ,'--','linewidth',lw,'color',color14); hold on
plot(P5.t/5,squeeze(P5.MoM.upup(i,j,:).^0.5),'-.','linewidth',lw,'color',color5); hold on


%%

close all
fig=figure('units','inches','outerposition',[1 1 12+1 6+1],'color','w');
iFig = fig.Number;
subplot(411)
plot(P1.t/1,abs(squeeze(M1.xpxp(i,j,:).^0.5)-squeeze(P1.MoM.xpxp(i,j,:).^0.5))./max(abs(squeeze(M1.xpxp(i,j,:).^0.5))),'-','linewidth',lw,'color',color1); hold on
plot(P2.t/2,abs(squeeze(M2.xpxp(i,j,:).^0.5)-squeeze(P2.MoM.xpxp(i,j,:).^0.5))./max(abs(squeeze(M2.xpxp(i,j,:).^0.5))),'-','linewidth',lw,'color',color2); hold on
plot(P3.t/3,abs(squeeze(M3.xpxp(i,j,:).^0.5)-squeeze(P3.MoM.xpxp(i,j,:).^0.5))./max(abs(squeeze(M3.xpxp(i,j,:).^0.5))),'-','linewidth',lw,'color',color3); hold on
plot(P4.t/4,abs(squeeze(M4.xpxp(i,j,:).^0.5)-squeeze(P4.MoM.xpxp(i,j,:).^0.5))./max(abs(squeeze(M4.xpxp(i,j,:).^0.5))),'-','linewidth',lw,'color',color4); hold on
plot(P5.t/5,abs(squeeze(M5.xpxp(i,j,:).^0.5)-squeeze(P5.MoM.xpxp(i,j,:).^0.5))./max(abs(squeeze(M5.xpxp(i,j,:).^0.5))),'-','linewidth',lw,'color',color5); hold on
leg = legend('$St=1$','$St=2$','$St=3$','$St=4$','$St=5$', ...
    'interpreter','latex','fontsize',0.7*fs,'location','northwest');
leg.Box = 'off';
leg.NumColumns = 5;
ylabel('$\varepsilon(\overline{Y}_1)$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ax.XAxis.TickLabels = [];
ylim([0 0.06])
subplot(412)
plot(P1.t/1,abs(squeeze(M1.mean_up(i,j,:).^0.5)-squeeze(P1.MoM.mean_up(i,j,:).^0.5))./max(abs(squeeze(M1.mean_up(i,j,:).^0.5))),'-','linewidth',lw,'color',color1); hold on
plot(P2.t/2,abs(squeeze(M2.mean_up(i,j,:).^0.5)-squeeze(P2.MoM.mean_up(i,j,:).^0.5))./max(abs(squeeze(M2.mean_up(i,j,:).^0.5))),'-','linewidth',lw,'color',color2); hold on
plot(P3.t/3,abs(squeeze(M3.mean_up(i,j,:).^0.5)-squeeze(P3.MoM.mean_up(i,j,:).^0.5))./max(abs(squeeze(M3.mean_up(i,j,:).^0.5))),'-','linewidth',lw,'color',color3); hold on
plot(P4.t/4,abs(squeeze(M4.mean_up(i,j,:).^0.5)-squeeze(P4.MoM.mean_up(i,j,:).^0.5))./max(abs(squeeze(M4.mean_up(i,j,:).^0.5))),'-','linewidth',lw,'color',color4); hold on
plot(P5.t/5,abs(squeeze(M5.mean_up(i,j,:).^0.5)-squeeze(P5.MoM.mean_up(i,j,:).^0.5))./max(abs(squeeze(M5.mean_up(i,j,:).^0.5))),'-','linewidth',lw,'color',color5); hold on
ylabel('$\varepsilon(\overline{V}_1)$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ax.XAxis.TickLabels = [];
ylim([0 0.02])
subplot(413)
plot(P1.t/1,abs(squeeze(M1.xpxp(i,j,:).^0.5)-squeeze(P1.MoM.xpxp(i,j,:).^0.5))./max(abs(squeeze(M1.xpxp(i,j,:).^0.5))),'-','linewidth',lw,'color',color1); hold on
plot(P2.t/2,abs(squeeze(M2.xpxp(i,j,:).^0.5)-squeeze(P2.MoM.xpxp(i,j,:).^0.5))./max(abs(squeeze(M2.xpxp(i,j,:).^0.5))),'-','linewidth',lw,'color',color2); hold on
plot(P3.t/3,abs(squeeze(M3.xpxp(i,j,:).^0.5)-squeeze(P3.MoM.xpxp(i,j,:).^0.5))./max(abs(squeeze(M3.xpxp(i,j,:).^0.5))),'-','linewidth',lw,'color',color3); hold on
plot(P4.t/4,abs(squeeze(M4.xpxp(i,j,:).^0.5)-squeeze(P4.MoM.xpxp(i,j,:).^0.5))./max(abs(squeeze(M4.xpxp(i,j,:).^0.5))),'-','linewidth',lw,'color',color4); hold on
plot(P5.t/5,abs(squeeze(M5.xpxp(i,j,:).^0.5)-squeeze(P5.MoM.xpxp(i,j,:).^0.5))./max(abs(squeeze(M5.xpxp(i,j,:).^0.5))),'-','linewidth',lw,'color',color5); hold on
ylabel('$\varepsilon(\sigma_{Y_1})$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ax.XAxis.TickLabels = [];
ylim([0 0.06])
subplot(414)
plot(P1.t/1,abs(squeeze(M1.upup(i,j,:).^0.5)-squeeze(P1.MoM.upup(i,j,:).^0.5))./max(abs(squeeze(M1.upup(i,j,:).^0.5))),'-','linewidth',lw,'color',color1); hold on
plot(P2.t/2,abs(squeeze(M2.upup(i,j,:).^0.5)-squeeze(P2.MoM.upup(i,j,:).^0.5))./max(abs(squeeze(M2.upup(i,j,:).^0.5))),'-','linewidth',lw,'color',color2); hold on
plot(P3.t/3,abs(squeeze(M3.upup(i,j,:).^0.5)-squeeze(P3.MoM.upup(i,j,:).^0.5))./max(abs(squeeze(M3.upup(i,j,:).^0.5))),'-','linewidth',lw,'color',color3); hold on
plot(P4.t/4,abs(squeeze(M4.upup(i,j,:).^0.5)-squeeze(P4.MoM.upup(i,j,:).^0.5))./max(abs(squeeze(M4.upup(i,j,:).^0.5))),'-','linewidth',lw,'color',color4); hold on
plot(P5.t/5,abs(squeeze(M5.upup(i,j,:).^0.5)-squeeze(P5.MoM.upup(i,j,:).^0.5))./max(abs(squeeze(M5.upup(i,j,:).^0.5))),'-','linewidth',lw,'color',color5); hold on
ylabel('$\varepsilon(\sigma_{V_1})$','interpreter','latex','fontsize',fs)
xlabel('$t/St$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylim([0 1])
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'sine_errors';
%     savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg


%%

close all
subplot(211)
plot(P1.t/1,abs(squeeze(M1.xpxp(i,j,:).^0.5)-squeeze(P1.MoM.xpxp(i,j,:).^0.5))./max(abs(squeeze(M1.xpxp(i,j,:).^0.5))),'-','linewidth',lw,'color',color1); hold on
plot(P2.t/2,abs(squeeze(M2.xpxp(i,j,:).^0.5)-squeeze(P2.MoM.xpxp(i,j,:).^0.5))./max(abs(squeeze(M2.xpxp(i,j,:).^0.5))),'-','linewidth',lw,'color',color2); hold on
plot(P3.t/3,abs(squeeze(M3.xpxp(i,j,:).^0.5)-squeeze(P3.MoM.xpxp(i,j,:).^0.5))./max(abs(squeeze(M3.xpxp(i,j,:).^0.5))),'-','linewidth',lw,'color',color3); hold on
plot(P4.t/4,abs(squeeze(M4.xpxp(i,j,:).^0.5)-squeeze(P4.MoM.xpxp(i,j,:).^0.5))./max(abs(squeeze(M4.xpxp(i,j,:).^0.5))),'-','linewidth',lw,'color',color4); hold on
plot(P5.t/5,abs(squeeze(M5.xpxp(i,j,:).^0.5)-squeeze(P5.MoM.xpxp(i,j,:).^0.5))./max(abs(squeeze(M5.xpxp(i,j,:).^0.5))),'-','linewidth',lw,'color',color5); hold on
subplot(212)
plot(P1.t/1,abs(squeeze(M1.mean_up(i,j,:).^0.5)-squeeze(P1.MoM.mean_up(i,j,:).^0.5))./max(abs(squeeze(M1.mean_up(i,j,:).^0.5))),'-','linewidth',lw,'color',color1); hold on
plot(P2.t/2,abs(squeeze(M2.mean_up(i,j,:).^0.5)-squeeze(P2.MoM.mean_up(i,j,:).^0.5))./max(abs(squeeze(M2.mean_up(i,j,:).^0.5))),'-','linewidth',lw,'color',color2); hold on
plot(P3.t/3,abs(squeeze(M3.mean_up(i,j,:).^0.5)-squeeze(P3.MoM.mean_up(i,j,:).^0.5))./max(abs(squeeze(M3.mean_up(i,j,:).^0.5))),'-','linewidth',lw,'color',color3); hold on
plot(P4.t/4,abs(squeeze(M4.mean_up(i,j,:).^0.5)-squeeze(P4.MoM.mean_up(i,j,:).^0.5))./max(abs(squeeze(M4.mean_up(i,j,:).^0.5))),'-','linewidth',lw,'color',color4); hold on
plot(P5.t/5,abs(squeeze(M5.mean_up(i,j,:).^0.5)-squeeze(P5.MoM.mean_up(i,j,:).^0.5))./max(abs(squeeze(M5.mean_up(i,j,:).^0.5))),'-','linewidth',lw,'color',color5); hold on




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
tLim       = [0 20];

% Forcing and initial condition
mean_a1  = 1;
sigma_a1 = 0.1;
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
dataf1.mean_a1  = mean_a1;
dataf1.a1       = a1;
dataf1.a1a1     = sigma_a1^2;
dataf1.a1a1a1   = cm3_a1;
dataf1.a1a1a1a1 = cm4_a1;
P1 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,[tLim(1) taup*tLim(2)],[taup*nt 2], ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

% -------------------------------------------------------------------------
taup = 3;
clear xp0 yp0 up0 vp0
xp0(1,1) = xp00;
yp0(1,1) = yp00;
up0(1,1) = up00;
vp0(1,1) = vp00;
dataf1.mean_a1  = mean_a1;
dataf1.a1       = a1;
dataf1.a1a1     = sigma_a1^2;
dataf1.a1a1a1   = cm3_a1;
dataf1.a1a1a1a1 = cm4_a1;
P3 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,[tLim(1) taup*tLim(2)],[taup*nt 2], ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

% -------------------------------------------------------------------------
taup = 5;
clear xp0 yp0 up0 vp0
xp0(1,1) = xp00;
yp0(1,1) = yp00;
up0(1,1) = up00;
vp0(1,1) = vp00;
dataf1.mean_a1  = mean_a1;
dataf1.a1       = a1;
dataf1.a1a1     = sigma_a1^2;
dataf1.a1a1a1   = cm3_a1;
dataf1.a1a1a1a1 = cm4_a1;
P5 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,[tLim(1) taup*tLim(2)],[taup*nt 2], ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup);


%% SPARSE-R ===============================================================

% -------------------------------------------------------------------------
taup     = 1;
mean_xp0 = xp00;
mean_up0 = up00;
mean_a   = mean_a1;
aa       = sigma_a1^2;
S1 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,[tLim(1) taup*tLim(2)],[taup*nt 2],flowType,dataFlow,timeMethod,1,taup);
S12 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,[tLim(1) taup*tLim(2)],[taup*nt 2],flowType,dataFlow,timeMethod,2,taup);
S13 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,[tLim(1) taup*tLim(2)],[taup*nt 2],flowType,dataFlow,timeMethod,3,taup);

% -------------------------------------------------------------------------
taup = 3;
S3 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,[tLim(1) taup*tLim(2)],[taup*nt 2],flowType,dataFlow,timeMethod,1,taup);
S32 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,[tLim(1) taup*tLim(2)],[taup*nt 2],flowType,dataFlow,timeMethod,2,taup);
S33 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,[tLim(1) taup*tLim(2)],[taup*nt 2],flowType,dataFlow,timeMethod,3,taup);

% -------------------------------------------------------------------------
taup = 5;
S5 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,[tLim(1) taup*tLim(2)],[taup*nt 2],flowType,dataFlow,timeMethod,1,taup);
S52 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,[tLim(1) taup*tLim(2)],[taup*nt 2],flowType,dataFlow,timeMethod,2,taup);
S53 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,[tLim(1) taup*tLim(2)],[taup*nt 2],flowType,dataFlow,timeMethod,3,taup);


%% FIGURE ERRORS MoP ======================================================

close all
i = 1;
j = 1;

e1_mean_xp  = abs(squeeze(P1.MoM.mean_xp(i,j,:))-squeeze(M1.mean_xp(i,j,:)))/norm(squeeze(M1.mean_xp(i,j,:)),Inf);
e3_mean_xp  = abs(squeeze(P3.MoM.mean_xp(i,j,:))-squeeze(M3.mean_xp(i,j,:)))/norm(squeeze(M3.mean_xp(i,j,:)),Inf);
e5_mean_xp  = abs(squeeze(P5.MoM.mean_xp(i,j,:))-squeeze(M5.mean_xp(i,j,:)))/norm(squeeze(M5.mean_xp(i,j,:)),Inf);
e1_mean_up  = abs(squeeze(P1.MoM.mean_up(i,j,:))-squeeze(M1.mean_up(i,j,:)))/norm(squeeze(M1.mean_up(i,j,:)),Inf);
e3_mean_up  = abs(squeeze(P3.MoM.mean_up(i,j,:))-squeeze(M3.mean_up(i,j,:)))/norm(squeeze(M3.mean_up(i,j,:)),Inf);
e5_mean_up  = abs(squeeze(P5.MoM.mean_up(i,j,:))-squeeze(M5.mean_up(i,j,:)))/norm(squeeze(M5.mean_up(i,j,:)),Inf);
e1_sigma_xp = abs(squeeze(P1.MoM.xpxp(i,j,:).^0.5)-squeeze(M1.xpxp(i,j,:).^0.5))/norm(squeeze(M1.xpxp(i,j,:).^0.5),Inf);
e3_sigma_xp = abs(squeeze(P3.MoM.xpxp(i,j,:).^0.5)-squeeze(M3.xpxp(i,j,:).^0.5))/norm(squeeze(M3.xpxp(i,j,:).^0.5),Inf);
e5_sigma_xp = abs(squeeze(P5.MoM.xpxp(i,j,:).^0.5)-squeeze(M5.xpxp(i,j,:).^0.5))/norm(squeeze(M5.xpxp(i,j,:).^0.5),Inf);
e1_xpup     = abs(squeeze(P1.MoM.xpup(i,j,:))-squeeze(M1.xpup(i,j,:)))/norm(squeeze(M1.xpup(i,j,:)),Inf);
e3_xpup     = abs(squeeze(P3.MoM.xpup(i,j,:))-squeeze(M3.xpup(i,j,:)))/norm(squeeze(M3.xpup(i,j,:)),Inf);
e5_xpup     = abs(squeeze(P5.MoM.xpup(i,j,:))-squeeze(M5.xpup(i,j,:)))/norm(squeeze(M5.xpup(i,j,:)),Inf);
e1_sigma_up = abs(squeeze(P1.MoM.upup(i,j,:).^0.5)-squeeze(M1.upup(i,j,:).^0.5))/norm(squeeze(M1.upup(i,j,:).^0.5),Inf);
e3_sigma_up = abs(squeeze(P3.MoM.upup(i,j,:).^0.5)-squeeze(M3.upup(i,j,:).^0.5))/norm(squeeze(M3.upup(i,j,:).^0.5),Inf);
e5_sigma_up = abs(squeeze(P5.MoM.upup(i,j,:).^0.5)-squeeze(M5.upup(i,j,:).^0.5))/norm(squeeze(M5.upup(i,j,:).^0.5),Inf);

fig      = figure('units','inches','outerposition',[1 1 12+1 5+1],'color','w');
iFig     = fig.Number;
[ha,pos] = tight_subplot(4,1,[0.065 0.01],0.15,0.08);

axes(ha(1)) % close all
% subplot(411)
plot(P1.t/1,e1_mean_xp,'-','linewidth',lw,'color',color1); hold on
plot(P3.t/3,e3_mean_xp,'-','linewidth',lw,'color',color2); hold on
plot(P5.t/5,e5_mean_xp,'-','linewidth',lw,'color',color3); hold on
leg = legend('$St=1$','$St=3$','$St=5$', ...
    'interpreter','latex','fontsize',0.7*fs,'location','northwest');
leg.Box = 'off';
leg.NumColumns = 5;
ylabel('$\varepsilon(\overline{Y}_1)$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ax.XAxis.TickLabels = [];
ax.YAxis.Exponent = -3;
ylim([0 0.002])
% subplot(412)
axes(ha(2)) % close all
plot(P1.t/1,e1_mean_up,'-','linewidth',lw,'color',color1); hold on
plot(P3.t/3,e3_mean_up,'-','linewidth',lw,'color',color2); hold on
plot(P5.t/5,e5_mean_up,'-','linewidth',lw,'color',color3); hold on
ylabel('$\varepsilon(\overline{V}_1)$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ax.XAxis.TickLabels = [];
ax.YAxis.Exponent = -2;
ylim([0 0.04])
% subplot(413)
axes(ha(3)) % close all
plot(P1.t/1,e1_sigma_xp,'-','linewidth',lw,'color',color1); hold on
plot(P3.t/3,e3_sigma_xp,'-','linewidth',lw,'color',color2); hold on
plot(P5.t/5,e5_sigma_xp,'-','linewidth',lw,'color',color3); hold on
ylabel('$\varepsilon(\sigma_{Y_1})$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ax.XAxis.TickLabels = [];
ax.YAxis.Exponent = -2;
ylim([0 0.04])
% subplot(414)
axes(ha(4)) % close all
plot(P1.t/1,e1_sigma_up,'-','linewidth',lw,'color',color1); hold on
plot(P3.t/3,e3_sigma_up,'-','linewidth',lw,'color',color2); hold on
plot(P5.t/5,e5_sigma_up,'-','linewidth',lw,'color',color3); hold on
ylabel('$\varepsilon(\sigma_{V_1})$','interpreter','latex','fontsize',fs)
xlabel('$t/St$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ax.YAxis.Exponent = 0;
ylim([0 1])
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'sine_errors';
    savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg


%% FIGURE ERRORS SPARSE-R M=1 =============================================

% close all

se1_mean_xp  = abs(S1.mean_xp-squeeze(M1.mean_xp(i,j,:)))/norm(squeeze(M1.mean_xp(i,j,:)),Inf);
se3_mean_xp  = abs(S3.mean_xp-squeeze(M3.mean_xp(i,j,:)))/norm(squeeze(M3.mean_xp(i,j,:)),Inf);
se5_mean_xp  = abs(S5.mean_xp-squeeze(M5.mean_xp(i,j,:)))/norm(squeeze(M5.mean_xp(i,j,:)),Inf);
se1_mean_up  = abs(S1.mean_up-squeeze(M1.mean_up(i,j,:)))/norm(squeeze(M1.mean_up(i,j,:)),Inf);
se3_mean_up  = abs(S3.mean_up-squeeze(M3.mean_up(i,j,:)))/norm(squeeze(M3.mean_up(i,j,:)),Inf);
se5_mean_up  = abs(S5.mean_up-squeeze(M5.mean_up(i,j,:)))/norm(squeeze(M5.mean_up(i,j,:)),Inf);
se1_sigma_xp = abs(S1.xpxp.^0.5-squeeze(M1.xpxp(i,j,:).^0.5))/norm(squeeze(M1.xpxp(i,j,:).^0.5),Inf);
se3_sigma_xp = abs(S3.xpxp.^0.5-squeeze(M3.xpxp(i,j,:).^0.5))/norm(squeeze(M3.xpxp(i,j,:).^0.5),Inf);
se5_sigma_xp = abs(S5.xpxp.^0.5-squeeze(M5.xpxp(i,j,:).^0.5))/norm(squeeze(M5.xpxp(i,j,:).^0.5),Inf);
se1_xpup     = abs(S1.xpup-squeeze(M1.xpup(i,j,:)))/norm(squeeze(M1.xpup(i,j,:)),Inf);
se3_xpup     = abs(S3.xpup-squeeze(M3.xpup(i,j,:)))/norm(squeeze(M3.xpup(i,j,:)),Inf);
se5_xpup     = abs(S5.xpup-squeeze(M5.xpup(i,j,:)))/norm(squeeze(M5.xpup(i,j,:)),Inf);
se1_sigma_up = abs(S1.upup.^0.5-squeeze(M1.upup(i,j,:).^0.5))/norm(squeeze(M1.upup(i,j,:).^0.5),Inf);
se3_sigma_up = abs(S3.upup.^0.5-squeeze(M3.upup(i,j,:).^0.5))/norm(squeeze(M3.upup(i,j,:).^0.5),Inf);
se5_sigma_up = abs(S5.upup.^0.5-squeeze(M5.upup(i,j,:).^0.5))/norm(squeeze(M5.upup(i,j,:).^0.5),Inf);

fig      = figure('units','inches','outerposition',[1 1 12+1 5+1],'color','w');
iFig     = fig.Number;
[ha,pos] = tight_subplot(4,1,[0.065 0.01],0.15,0.08);

axes(ha(1)) % close all
% subplot(411)
plot(P1.t/1,se1_mean_xp,'-','linewidth',lw,'color',color1); hold on
plot(P3.t/3,se3_mean_xp,'-','linewidth',lw,'color',color2); hold on
plot(P5.t/5,se5_mean_xp,'-','linewidth',lw,'color',color3); hold on
leg = legend('$St=1$','$St=3$','$St=5$', ...
    'interpreter','latex','fontsize',0.7*fs,'location','northwest');
leg.Box = 'off';
leg.NumColumns = 5;
ylabel('$\varepsilon(\overline{Y}_1)$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ax.XAxis.TickLabels = [];
ax.YAxis.Exponent = -4;
ylim([0 0.0002])
% subplot(412)
axes(ha(2)) % close all
plot(P1.t/1,se1_mean_up,'-','linewidth',lw,'color',color1); hold on
plot(P3.t/3,se3_mean_up,'-','linewidth',lw,'color',color2); hold on
plot(P5.t/5,se5_mean_up,'-','linewidth',lw,'color',color3); hold on
ylabel('$\varepsilon(\overline{V}_1)$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ax.XAxis.TickLabels = [];
ax.YAxis.Exponent = -2;
ylim([0 0.02])
% subplot(413)
axes(ha(3)) % close all
plot(P1.t/1,se1_sigma_xp,'-','linewidth',lw,'color',color1); hold on
plot(P3.t/3,se3_sigma_xp,'-','linewidth',lw,'color',color2); hold on
plot(P5.t/5,se5_sigma_xp,'-','linewidth',lw,'color',color3); hold on
ylabel('$\varepsilon(\sigma_{Y_1})$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ax.XAxis.TickLabels = [];
ax.YAxis.Exponent = -1;
ax.YAxis.TickValues = linspace(0,0.6,3);
ylim([0 0.6])
% subplot(414)
axes(ha(4)) % close all
plot(P1.t/1,se1_sigma_up,'-','linewidth',lw,'color',color1); hold on
plot(P3.t/3,se3_sigma_up,'-','linewidth',lw,'color',color2); hold on
plot(P5.t/5,se5_sigma_up,'-','linewidth',lw,'color',color3); hold on
ylabel('$\varepsilon(\sigma_{V_1})$','interpreter','latex','fontsize',fs)
xlabel('$t/St$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ax.YAxis.Exponent = -1;
ylim([0 0.4])
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'sine_errorsSPARSER';
    savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg


%% FIGURE ERRORS SPARSE-R M=3 =============================================

% close all

se1_mean_xp  = abs(S13.mean_xp-squeeze(M1.mean_xp(i,j,:)))/norm(squeeze(M1.mean_xp(i,j,:)),Inf);
se3_mean_xp  = abs(S33.mean_xp-squeeze(M3.mean_xp(i,j,:)))/norm(squeeze(M3.mean_xp(i,j,:)),Inf);
se5_mean_xp  = abs(S53.mean_xp-squeeze(M5.mean_xp(i,j,:)))/norm(squeeze(M5.mean_xp(i,j,:)),Inf);
se1_mean_up  = abs(S13.mean_up-squeeze(M1.mean_up(i,j,:)))/norm(squeeze(M1.mean_up(i,j,:)),Inf);
se3_mean_up  = abs(S33.mean_up-squeeze(M3.mean_up(i,j,:)))/norm(squeeze(M3.mean_up(i,j,:)),Inf);
se5_mean_up  = abs(S53.mean_up-squeeze(M5.mean_up(i,j,:)))/norm(squeeze(M5.mean_up(i,j,:)),Inf);
se1_sigma_xp = abs(S13.xpxp.^0.5-squeeze(M1.xpxp(i,j,:).^0.5))/norm(squeeze(M1.xpxp(i,j,:).^0.5),Inf);
se3_sigma_xp = abs(S33.xpxp.^0.5-squeeze(M3.xpxp(i,j,:).^0.5))/norm(squeeze(M3.xpxp(i,j,:).^0.5),Inf);
se5_sigma_xp = abs(S53.xpxp.^0.5-squeeze(M5.xpxp(i,j,:).^0.5))/norm(squeeze(M5.xpxp(i,j,:).^0.5),Inf);
se1_xpup     = abs(S13.xpup-squeeze(M1.xpup(i,j,:)))/norm(squeeze(M1.xpup(i,j,:)),Inf);
se3_xpup     = abs(S33.xpup-squeeze(M3.xpup(i,j,:)))/norm(squeeze(M3.xpup(i,j,:)),Inf);
se5_xpup     = abs(S53.xpup-squeeze(M5.xpup(i,j,:)))/norm(squeeze(M5.xpup(i,j,:)),Inf);
se1_sigma_up = abs(S13.upup.^0.5-squeeze(M1.upup(i,j,:).^0.5))/norm(squeeze(M1.upup(i,j,:).^0.5),Inf);
se3_sigma_up = abs(S33.upup.^0.5-squeeze(M3.upup(i,j,:).^0.5))/norm(squeeze(M3.upup(i,j,:).^0.5),Inf);
se5_sigma_up = abs(S53.upup.^0.5-squeeze(M5.upup(i,j,:).^0.5))/norm(squeeze(M5.upup(i,j,:).^0.5),Inf);

fig      = figure('units','inches','outerposition',[1 1 12+1 5+1],'color','w');
iFig     = fig.Number;
[ha,pos] = tight_subplot(4,1,[0.065 0.01],0.15,0.08);

axes(ha(1)) % close all
% subplot(411)
plot(P1.t/1,se1_mean_xp,'-','linewidth',lw,'color',color1); hold on
plot(P3.t/3,se3_mean_xp,'-','linewidth',lw,'color',color2); hold on
plot(P5.t/5,se5_mean_xp,'-','linewidth',lw,'color',color3); hold on
leg = legend('$St=1$','$St=3$','$St=5$', ...
    'interpreter','latex','fontsize',0.7*fs,'location','northwest');
leg.Box = 'off';
leg.NumColumns = 5;
ylabel('$\varepsilon(\overline{Y}_1)$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ax.XAxis.TickLabels = [];
% ylim([0 0.002])
% subplot(412)
axes(ha(2)) % close all
plot(P1.t/1,se1_mean_up,'-','linewidth',lw,'color',color1); hold on
plot(P3.t/3,se3_mean_up,'-','linewidth',lw,'color',color2); hold on
plot(P5.t/5,se5_mean_up,'-','linewidth',lw,'color',color3); hold on
ylabel('$\varepsilon(\overline{V}_1)$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ax.XAxis.TickLabels = [];
% ylim([0 0.04])
axes(ha(3)) % close all
% subplot(413)
plot(P1.t/1,se1_sigma_xp,'-','linewidth',lw,'color',color1); hold on
plot(P3.t/3,se3_sigma_xp,'-','linewidth',lw,'color',color2); hold on
plot(P5.t/5,se5_sigma_xp,'-','linewidth',lw,'color',color3); hold on
ylabel('$\varepsilon(\sigma_{Y_1})$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ax.XAxis.TickLabels = [];
% ylim([0 0.04])
% subplot(414)
axes(ha(4)) % close all
plot(P1.t/1,se1_sigma_up,'-','linewidth',lw,'color',color1); hold on
plot(P3.t/3,se3_sigma_up,'-','linewidth',lw,'color',color2); hold on
plot(P5.t/5,se5_sigma_up,'-','linewidth',lw,'color',color3); hold on
ylabel('$\varepsilon(\sigma_{V_1})$','interpreter','latex','fontsize',fs)
xlabel('$t/St$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
% ylim([0 1])
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'sine_errorsSPARSER3';
    savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg




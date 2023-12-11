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

% -------------------------------------------------------------------------
taup = 5;
ns   = 10000; % Numer of samples per cloud in MC-PSIC
a1 = linspace(mean_a1-0.5*sqrt(12)*sigma_a1,mean_a1+0.5*sqrt(12)*sigma_a1,ns)';
for k=1:ns
    xp0(1,1,k)       = xp00;
    yp0(1,1,k)       = yp00;
    up0(1,1,k)       = up00;
    vp0(1,1,k)       = vp00;
    dataf1.a1(1,1,k) = a1(k);
end
M5 = FUNC_solverMCPSIC2D_inertial(xp0,yp0,up0,vp0,tLim,nt,flowType,dataFlow,f1Type,dataf1,timeMethod,taup);


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
P1 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,tLim,nt, ...
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
P5 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,tLim,nt, ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup);


%% SPARSER ================================================================

% -------------------------------------------------------------------------
taup     = 1;
mean_xp0 = xp00;
mean_up0 = up00;
mean_a   = mean_a1;
aa       = sigma_a1^2;
S11 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,1,taup);
S12 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,2,taup);
S13 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,3,taup);
S14 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,4,taup);
S15 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,5,taup);

% -------------------------------------------------------------------------
taup     = 5;
mean_xp0 = xp00;
mean_up0 = up00;
mean_a   = mean_a1;
aa       = sigma_a1^2;
S51 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,1,taup);
S52 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,2,taup);
S53 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,3,taup);
S54 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,4,taup);
S55 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,5,taup);




close all

subplot(221)
plot(M1.t,squeeze(M1.mean_xp(1,1,:))); hold on
plot(P1.t,squeeze(P1.MoM.mean_xp(1,1,:)),'--'); hold on
plot(S11.t,S11.mean_xp,'-.')
plot(M1.t,squeeze(M1.mean_up(1,1,:))); hold on
plot(P1.t,squeeze(P1.MoM.mean_up(1,1,:)),'--'); hold on
plot(S11.t,S11.mean_up,'-.')

subplot(222)
plot(M1.t,squeeze(M1.xpxp(1,1,:)).^0.5); hold on
plot(P1.t,squeeze(P1.MoM.xpxp(1,1,:)).^0.5,'--'); hold on
plot(S11.t,S11.xpxp.^0.5,'-.')
plot(M1.t,squeeze(M1.upup(1,1,:)).^0.5); hold on
plot(P1.t,squeeze(P1.MoM.upup(1,1,:)).^0.5,'--'); hold on
plot(S11.t,S11.upup.^0.5,'-.')

subplot(223)
plot(M5.t,squeeze(M5.mean_xp(1,1,:))); hold on
plot(P5.t,squeeze(P5.MoM.mean_xp(1,1,:)),'--'); hold on
plot(S51.t,S51.mean_xp,'-.')
plot(M5.t,squeeze(M5.mean_up(1,1,:))); hold on
plot(P5.t,squeeze(P5.MoM.mean_up(1,1,:)),'--'); hold on
plot(S51.t,S51.mean_up,'-.')

subplot(224)
plot(M5.t,squeeze(M5.xpxp(1,1,:)).^0.5); hold on
plot(P5.t,squeeze(P5.MoM.xpxp(1,1,:)).^0.5,'--'); hold on
plot(S51.t,S51.xpxp.^0.5,'-.')
plot(M5.t,squeeze(M5.upup(1,1,:)).^0.5); hold on
plot(P5.t,squeeze(P5.MoM.upup(1,1,:)).^0.5,'--'); hold on
plot(S15.t,S51.upup.^0.5,'-.')


% figure % close all
% plot([1 2 3 4 5],[e_xp_S1 e_xp_S2 e_xp_S3 e_xp_S4 e_xp_S5],'o--'); hold on
% plot([1 5],[e_xp_P1 e_xp_P1])


%% PLOTS ==================================================================

close all
i  = 1;
j  = 1;
nn = round(nt(1)/3);

% -------------------------------------------------------------------------
figure(iFig) % close all
p1(1) = plot(P1.t,squeeze(P1.o0.xp(i,j,:)),'linewidth',lw,'color',color14); hold on
p1(2) = plot(P1.t,squeeze(P1.o1.xp(i,j,:)),'--','linewidth',lw,'color',color14); hold on
p1(3) = plot(P1.t,squeeze(P1.o2.xp(i,j,:)),'-.','linewidth',lw,'color',color14); hold on
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylim([-4 12])
yyaxis right
p1(4) = plot(P1.t,squeeze(P1.o0.up(i,j,:)),'linewidth',lw,'color',color8); hold on
p1(5) = plot(P1.t,squeeze(P1.o1.up(i,j,:)),'--','linewidth',lw,'color',color8); hold on
p1(6) = plot(P1.t,squeeze(P1.o2.up(i,j,:)),'-.','linewidth',lw,'color',color8); hold on
ylim([-2 3])
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylabel('$v_{1(0)}, \ v_{1(2)}, \ v_{1(2)}$','interpreter','latex','fontsize',fs)
xlabel('$t$','interpreter','latex','fontsize',fs)
yyaxis left
ylabel('$y_{1(0)}, \ y_{1(2)}, \ y_{1(2)}$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
leg1{1} = '$y_{1(0)}$';
leg1{2} = '$y_{1(1)}$';
leg1{3} = '$y_{1(2)}$';
leg1{4} = '$v_{1(0)}$';
leg1{5} = '$v_{1(1)}$';
leg1{6} = '$v_{1(2)}$';
leg     = legend(p1,leg1,'interpreter','latex','fontsize',0.8*fs,'location','north');
leg.Box = 'off';
leg.NumColumns = 2;
ax.XTick = 0:2:10;
ax.YAxis(1).TickValues = -20:2:20;
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'sine_orders_St1';
%     savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg


% -------------------------------------------------------------------------
figure(iFig) % close all
p1(1) = plot(P5.t,squeeze(P5.o0.xp(i,j,:)),'linewidth',lw,'color',color14); hold on
p1(2) = plot(P5.t,squeeze(P5.o1.xp(i,j,:)),'--','linewidth',lw,'color',color14); hold on
p1(3) = plot(P5.t,squeeze(P5.o2.xp(i,j,:)),'-.','linewidth',lw,'color',color14); hold on
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylim([-4 12])
yyaxis right
p1(4) = plot(P5.t,squeeze(P5.o0.up(i,j,:)),'linewidth',lw,'color',color8); hold on
p1(5) = plot(P5.t,squeeze(P5.o1.up(i,j,:)),'--','linewidth',lw,'color',color8); hold on
p1(6) = plot(P5.t,squeeze(P5.o2.up(i,j,:)),'-.','linewidth',lw,'color',color8); hold on
ylim([-2 3])
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylabel('$v_{1(0)}, \ v_{1(2)}, \ v_{1(2)}$','interpreter','latex','fontsize',fs)
xlabel('$t$','interpreter','latex','fontsize',fs)
yyaxis left
ylabel('$y_{1(0)}, \ y_{1(2)}, \ y_{1(2)}$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
leg1{1} = '$y_{1(0)}$';
leg1{2} = '$y_{1(1)}$';
leg1{3} = '$y_{1(2)}$';
leg1{4} = '$v_{1(0)}$';
leg1{5} = '$v_{1(1)}$';
leg1{6} = '$v_{1(2)}$';
% leg = legend(p1,leg1,'interpreter','latex','fontsize',0.8*fs,'location','north');
% leg.Box = 'off';
% leg.NumColumns = 2;
ax.XTick = 0:2:10;
ax.YAxis(1).TickValues = -20:2:20;
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'sine_orders_St5';
%     savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg


% -------------------------------------------------------------------------
figure(iFig) % close all
p1(1) = plot(M1.t,squeeze(M1.mean_xp(i,j,:)),'color',color9,'linewidth',lw); hold on
plot(M1.t,squeeze(M1.mean_xp(i,j,:))+2*squeeze(M1.xpxp(i,j,:).^0.5),'--','color',color9,'linewidth',lw); hold on
plot(M1.t,squeeze(M1.mean_xp(i,j,:))-2*squeeze(M1.xpxp(i,j,:).^0.5),'--','color',color9,'linewidth',lw); hold on
x1 = M1.t;
x2 = M1.t;                 % #create continuous x value array for plotting
y1 = squeeze(M1.mean_xp(i,j,:)) + 2*squeeze(M1.xpxp(i,j,:).^0.5);
y2 = squeeze(M1.mean_xp(i,j,:)) - 2*squeeze(M1.xpxp(i,j,:).^0.5);
X = [x1',fliplr(x2')];                 % #create continuous x value array for plotting
Y = [y1',fliplr(y2')];              %#create y values for out and then back
fill(X,Y,color9,'facealpha',0.2,'edgecolor','none');
ylim([0 10])
yyaxis right
p1(2) = plot(M1.t,squeeze(M1.mean_up(i,j,:)),'color',color8,'linewidth',lw); hold on
plot(M1.t,squeeze(M1.mean_up(i,j,:))+2*squeeze(M1.upup(i,j,:).^0.5),'--','color',color8,'linewidth',lw); hold on
plot(M1.t,squeeze(M1.mean_up(i,j,:))-2*squeeze(M1.upup(i,j,:).^0.5),'--','color',color8,'linewidth',lw); hold on
x1 = M1.t;
x2 = M1.t;                 % #create continuous x value array for plotting
y1 = squeeze(M1.mean_up(i,j,:)) + 2*squeeze(M1.upup(i,j,:).^0.5);
y2 = squeeze(M1.mean_up(i,j,:)) - 2*squeeze(M1.upup(i,j,:).^0.5);
X = [x1',fliplr(x2')];                 % #create continuous x value array for plotting
Y = [y1',fliplr(y2')];              %#create y values for out and then back
fill(X,Y,color8,'facealpha',0.2,'edgecolor','none'); 
ylim([0 1.5])
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylabel('$\overline{V}_1 \pm 2\sigma_{V_1}$','interpreter','latex','fontsize',fs)
xlabel('$t$','interpreter','latex','fontsize',fs)
ax.YTick = 0:0.5:1.5;
yyaxis left
ylabel('$\overline{Y}_1 \pm 2\sigma_{Y_1}$','interpreter','latex','fontsize',fs)
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
leg1{1} = '$\overline{Y}_1$';
leg1{2} = '$\overline{V}_2$';
leg = legend(p1,leg1,'interpreter','latex','fontsize',0.8*fs,'location','north');
leg.Box = 'off';
leg.NumColumns = 2;
ax.XTick = 0:2:10;
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'sine_bound_St1';
%     savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg

% -------------------------------------------------------------------------
figure(iFig) % close all
p1(1) = plot(M5.t,squeeze(M5.mean_xp(i,j,:)),'color',color9,'linewidth',lw); hold on
plot(M5.t,squeeze(M5.mean_xp(i,j,:))+2*squeeze(M5.xpxp(i,j,:).^0.5),'--','color',color9,'linewidth',lw); hold on
plot(M5.t,squeeze(M5.mean_xp(i,j,:))-2*squeeze(M5.xpxp(i,j,:).^0.5),'--','color',color9,'linewidth',lw); hold on
x1 = M5.t;
x2 = M5.t;                 % #create continuous x value array for plotting
y1 = squeeze(M5.mean_xp(i,j,:)) + 2*squeeze(M5.xpxp(i,j,:).^0.5);
y2 = squeeze(M5.mean_xp(i,j,:)) - 2*squeeze(M5.xpxp(i,j,:).^0.5);
X = [x1',fliplr(x2')];                 % #create continuous x value array for plotting
Y = [y1',fliplr(y2')];              %#create y values for out and then back
fill(X,Y,color9,'facealpha',0.2,'edgecolor','none');
ylim([0 10])
yyaxis right
p1(2) = plot(M5.t,squeeze(M5.mean_up(i,j,:)),'color',color8,'linewidth',lw); hold on
plot(M5.t,squeeze(M5.mean_up(i,j,:))+2*squeeze(M5.upup(i,j,:).^0.5),'--','color',color8,'linewidth',lw); hold on
plot(M5.t,squeeze(M5.mean_up(i,j,:))-2*squeeze(M5.upup(i,j,:).^0.5),'--','color',color8,'linewidth',lw); hold on
x1 = M5.t;
x2 = M5.t;                 % #create continuous x value array for plotting
y1 = squeeze(M5.mean_up(i,j,:)) + 2*squeeze(M5.upup(i,j,:).^0.5);
y2 = squeeze(M5.mean_up(i,j,:)) - 2*squeeze(M5.upup(i,j,:).^0.5);
X = [x1',fliplr(x2')];                 % #create continuous x value array for plotting
Y = [y1',fliplr(y2')];              %#create y values for out and then back
fill(X,Y,color8,'facealpha',0.2,'edgecolor','none'); 
ylim([0 1.5])
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylabel('$\overline{V}_1 \pm 2\sigma_{V_1}$','interpreter','latex','fontsize',fs)
xlabel('$t$','interpreter','latex','fontsize',fs)
ax.YTick = 0:0.5:1.5;
yyaxis left
ylabel('$\overline{Y}_1 \pm 2\sigma_{Y_1}$','interpreter','latex','fontsize',fs)
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
leg1{1} = '$\overline{Y}_1$';
leg1{2} = '$\overline{V}_2$';
% leg = legend(p1,leg1,'interpreter','latex','fontsize',0.8*fs,'location','north');
% leg.Box = 'off';
% leg.NumColumns = 2;
ax.XTick = 0:2:10;
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'sine_bound_St5';
%     savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg


% -------------------------------------------------------------------------
figure(iFig) % close all
p1(1) = plot(M1.t,squeeze(M1.xpxp(i,j,:).^0.5),'s-','linewidth',lw,'color',color14,'MarkerIndices',round(0.1*nn):nn:nt(1),'MarkerFaceColor',color14); hold on
p1(3) = plot(S11.t,squeeze(S11.xpxp.^0.5),'>--','linewidth',lw,'color',color8,'MarkerIndices',round(0.4*nn):nn:nt(1),'MarkerFaceColor',color8); hold on
p1(2) = plot(P1.t,squeeze(P1.MoM.xpxp(i,j,:).^0.5),'o-.','linewidth',lw,'color',color10,'MarkerIndices',round(0.7*nn):nn:nt(1),'MarkerFaceColor',color10); hold on
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylim([0 0.3])
yyaxis right
p1(4) = plot(M1.t,squeeze(M1.upup(i,j,:).^0.5),'d-','linewidth',lw,'color',color15,'MarkerIndices',round(0.1*nn):nn:nt(1),'MarkerFaceColor',color15); hold on
p1(6) = plot(S11.t,squeeze(S11.upup.^0.5),'^--','linewidth',lw,'color',color3,'MarkerIndices',round(0.4*nn):nn:nt(1),'MarkerFaceColor',color3); hold on
p1(5) = plot(P1.t,squeeze(P1.MoM.upup(i,j,:).^0.5),'*-.','linewidth',lw,'color',color9,'MarkerIndices',round(0.7*nn):nn:nt(1),'MarkerFaceColor',color9); hold on
ylim([0 0.1])
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylabel('$\sigma_{Y_1}$','interpreter','latex','fontsize',fs)
xlabel('$t$','interpreter','latex','fontsize',fs)
yyaxis left
ylabel('$\sigma_{V_1}$','interpreter','latex','fontsize',fs)
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
leg1{1} = '$\sigma_{Y_1}$, MC-PSIC';
leg1{2} = '$\sigma_{Y_1}$, MoM-MoP';
leg1{3} = '$\sigma_{Y_1}$, SPARSE-R';
leg1{4} = '$\sigma_{V_1}$, MC-PSIC';
leg1{5} = '$\sigma_{V_1}$, MoM-MoP';
leg1{6} = '$\sigma_{V_1}$, SPARSE-R';
leg = legend(p1,leg1,'interpreter','latex','fontsize',0.7*fs,'location','northeast');
leg.Box = 'off';
leg.NumColumns = 2;
ax.XTick = 0:2:10;
% ax.YTick = 0:0.02:0.12;
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'sine_cm2_St1';
%     savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg


% -------------------------------------------------------------------------
figure(iFig) % close all
p1(1) = plot(M5.t,squeeze(M5.xpxp(i,j,:).^0.5),'s-','linewidth',lw,'color',color14,'MarkerIndices',round(0.1*nn):nn:nt(1),'MarkerFaceColor',color14); hold on
p1(3) = plot(S51.t,squeeze(S51.xpxp.^0.5),'>--','linewidth',lw,'color',color8,'MarkerIndices',round(0.4*nn):nn:nt(1),'MarkerFaceColor',color8); hold on
p1(2) = plot(P5.t,squeeze(P5.MoM.xpxp(i,j,:).^0.5),'o-.','linewidth',lw,'color',color10,'MarkerIndices',round(0.7*nn):nn:nt(1),'MarkerFaceColor',color10); hold on
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylim([0 0.3])
yyaxis right
p1(4) = plot(M5.t,squeeze(M5.upup(i,j,:).^0.5),'d-','linewidth',lw,'color',color15,'MarkerIndices',round(0.1*nn):nn:nt(1),'MarkerFaceColor',color15); hold on
p1(6) = plot(S51.t,squeeze(S51.upup.^0.5),'^--','linewidth',lw,'color',color3,'MarkerIndices',round(0.4*nn):nn:nt(1),'MarkerFaceColor',color3); hold on
p1(5) = plot(P5.t,squeeze(P5.MoM.upup(i,j,:).^0.5),'*-.','linewidth',lw,'color',color9,'MarkerIndices',round(0.7*nn):nn:nt(1),'MarkerFaceColor',color9); hold on
ylim([0 0.1])
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylabel('$\sigma_{Y_1}$','interpreter','latex','fontsize',fs)
xlabel('$t$','interpreter','latex','fontsize',fs)
yyaxis left
ylabel('$\sigma_{V_1}$','interpreter','latex','fontsize',fs)
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
leg1{1} = '$\sigma_{Y_1}$, MC-PSIC';
leg1{2} = '$\sigma_{Y_1}$, MoM-MoP';
leg1{3} = '$\sigma_{Y_1}$, SPARSE-R';
leg1{4} = '$\sigma_{V_1}$, MC-PSIC';
leg1{5} = '$\sigma_{V_1}$, MoM-MoP';
leg1{6} = '$\sigma_{V_1}$, SPARSE-R';
% leg = legend(p1,leg1,'interpreter','latex','fontsize',0.8*fs,'location','northeast');
% leg.Box = 'off';
% leg.NumColumns = 1;
ax.XTick = 0:2:10;
% ax.YTick = 0:0.02:0.12;
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'sine_cm2_St5';
%     savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg





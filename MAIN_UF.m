%% ========================================================================
close all
clear
clc
setup
label = 'UF';
iFig  = 1;
saveFigs = 1;

% PARAMETERS --------------------------------------------------------------
timeMethod = 3;
flowType   = 0;
dataFlow   = [];
f1Type     = 3; % Stokes random
dataf1     = [];
nt         = [101 10];
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
a1   = linspace(mean_a1-0.5*sqrt(12)*sigma_a1,mean_a1+0.5*sqrt(12)*sigma_a1,ns)';
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
M5 = FUNC_solverMCPSIC2D_inertial(xp0,yp0,up0,vp0,tLim,nt,flowType,dataFlow,f1Type,dataf1,timeMethod,taup);


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
P1 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,tLim,nt, ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

% -------------------------------------------------------------------------
taup = 5;
P5 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,tLim,nt, ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup);


%% PLOTS St1 ==============================================================

close all
i = 1;
j = 1;
nn = 51;

% -------------------------------------------------------------------------
figure(iFig) % close all
p1(1) = plot(P1.t,squeeze(P1.o0.xp(i,j,:)),'linewidth',lw,'color',color14); hold on
p1(2) = plot(P1.t,squeeze(P1.o1.xp(i,j,:)),'--','linewidth',lw,'color',color14); hold on
p1(3) = plot(P1.t,squeeze(P1.o2.xp(i,j,:)),'-.','linewidth',lw,'color',color14); hold on
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylim([-2 10])
yyaxis right
p1(4) = plot(P1.t,squeeze(P1.o0.up(i,j,:)),'linewidth',lw,'color',color8); hold on
p1(5) = plot(P1.t,squeeze(P1.o1.up(i,j,:)),'--','linewidth',lw,'color',color8); hold on
p1(6) = plot(P1.t,squeeze(P1.o2.up(i,j,:)),'-.','linewidth',lw,'color',color8); hold on
ylim([-0.5 1])
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
leg = legend(p1,leg1,'interpreter','latex','fontsize',0.8*fs,'location','north');
leg.Box = 'off';
leg.NumColumns = 2;
ax.XTick = 0:2:10;
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'UF_orders_St1';
    savefig(iFig,thisName)
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
ylim([0 1.2])
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylabel('$\overline{V}_1 \pm 2\sigma_{V_1}$','interpreter','latex','fontsize',fs)
xlabel('$t$','interpreter','latex','fontsize',fs)
ax.YTick = 0:0.2:1.2;
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
    thisName = 'UF_bound_St1';
    savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg


% -------------------------------------------------------------------------
figure(iFig) % close all
p1(1) = plot(M1.t,squeeze(M1.xpxp(i,j,:).^0.5),'s-','linewidth',lw,'color',color14,'MarkerIndices',round(0.1*nn):nn:nt(1),'MarkerFaceColor',color14); hold on
p1(2) = plot(P1.t,squeeze(P1.xpxp(i,j,:).^0.5),'>--','linewidth',lw,'color',color8,'MarkerIndices',round(0.4*nn):nn:nt(1),'MarkerFaceColor',color8); hold on
p1(3) = plot(P1.t,squeeze(P1.MoM.xpxp(i,j,:).^0.5),'o-.','linewidth',lw,'color',color10,'MarkerIndices',round(0.7*nn):nn:nt(1),'MarkerFaceColor',color10); hold on
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylim([0 0.3])
yyaxis right
p1(4) = plot(M1.t,squeeze(M1.upup(i,j,:).^0.5),'d-','linewidth',lw,'color',color15,'MarkerIndices',round(0.1*nn):nn:nt(1),'MarkerFaceColor',color15); hold on
p1(5) = plot(P1.t,squeeze(P1.upup(i,j,:).^0.5),'^--','linewidth',lw,'color',color3,'MarkerIndices',round(0.4*nn):nn:nt(1),'MarkerFaceColor',color3); hold on
p1(6) = plot(P1.t,squeeze(P1.MoM.upup(i,j,:).^0.5),'*-.','linewidth',lw,'color',color9,'MarkerIndices',round(0.7*nn):nn:nt(1),'MarkerFaceColor',color9); hold on
ylim([0 0.04])
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
leg1{2} = '$\sigma_{Y_1}$, MC-MoP';
leg1{3} = '$\sigma_{Y_1}$, MoM-MoP';
leg1{4} = '$\sigma_{V_1}$, MC-PSIC';
leg1{5} = '$\sigma_{V_1}$, MC-MoP';
leg1{6} = '$\sigma_{V_1}$, MoM-MoP';
leg = legend(p1,leg1,'interpreter','latex','fontsize',0.8*fs,'location','northeast');
leg.Box = 'off';
leg.NumColumns = 1;
ax.XTick = 0:2:10;
% ax.YTick = 0:0.02:0.12;
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'UF_cm2_St1';
    savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg


%% St5 ====================================================================

close all
i = 1;
j = 1;
nn = 51;

% -------------------------------------------------------------------------
figure(iFig) % close all
p1(1) = plot(P5.t,squeeze(P5.o0.xp(i,j,:)),'linewidth',lw,'color',color14); hold on
p1(2) = plot(P5.t,squeeze(P5.o1.xp(i,j,:)),'--','linewidth',lw,'color',color14); hold on
p1(3) = plot(P5.t,squeeze(P5.o2.xp(i,j,:)),'-.','linewidth',lw,'color',color14); hold on
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylim([-2 10])
yyaxis right
p1(4) = plot(P5.t,squeeze(P5.o0.up(i,j,:)),'linewidth',lw,'color',color8); hold on
p1(5) = plot(P5.t,squeeze(P5.o1.up(i,j,:)),'--','linewidth',lw,'color',color8); hold on
p1(6) = plot(P5.t,squeeze(P5.o2.up(i,j,:)),'-.','linewidth',lw,'color',color8); hold on
ylim([-0.5 1])
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
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'UF_orders_St5';
    savefig(iFig,thisName)
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
ylim([0 1.2])
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylabel('$\overline{V}_1 \pm 2\sigma_{V_1}$','interpreter','latex','fontsize',fs)
xlabel('$t$','interpreter','latex','fontsize',fs)
ax.YTick = 0:0.2:1.2;
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
    thisName = 'UF_bound_St5';
    savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg


% -------------------------------------------------------------------------
figure(iFig) % close all
p1(1) = plot(M5.t,squeeze(M5.xpxp(i,j,:).^0.5),'s-','linewidth',lw,'color',color14,'MarkerIndices',round(0.1*nn):nn:nt(1),'MarkerFaceColor',color14); hold on
p1(2) = plot(P5.t,squeeze(P5.xpxp(i,j,:).^0.5),'>--','linewidth',lw,'color',color8,'MarkerIndices',round(0.4*nn):nn:nt(1),'MarkerFaceColor',color8); hold on
p1(3) = plot(P5.t,squeeze(P5.MoM.xpxp(i,j,:).^0.5),'o-.','linewidth',lw,'color',color10,'MarkerIndices',round(0.7*nn):nn:nt(1),'MarkerFaceColor',color10); hold on
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylim([0 0.3])
yyaxis right
p1(4) = plot(M5.t,squeeze(M5.upup(i,j,:).^0.5),'d-','linewidth',lw,'color',color15,'MarkerIndices',round(0.1*nn):nn:nt(1),'MarkerFaceColor',color15); hold on
p1(5) = plot(P5.t,squeeze(P5.upup(i,j,:).^0.5),'^--','linewidth',lw,'color',color3,'MarkerIndices',round(0.4*nn):nn:nt(1),'MarkerFaceColor',color3); hold on
p1(6) = plot(P5.t,squeeze(P5.MoM.upup(i,j,:).^0.5),'*-.','linewidth',lw,'color',color9,'MarkerIndices',round(0.7*nn):nn:nt(1),'MarkerFaceColor',color9); hold on
ylim([0 0.04])
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
leg1{2} = '$\sigma_{Y_1}$, MC-MoP';
leg1{3} = '$\sigma_{Y_1}$, MoM-MoP';
leg1{4} = '$\sigma_{V_1}$, MC-PSIC';
leg1{5} = '$\sigma_{V_1}$, MC-MoP';
leg1{6} = '$\sigma_{V_1}$, MoM-MoP';
% leg = legend(p1,leg1,'interpreter','latex','fontsize',0.8*fs,'location','east');
% leg.Box = 'off';
% leg.NumColumns = 1;
ax.XTick = 0:2:10;
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'UF_cm2_St5';
    savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg


%% CHECKING THE ANALYTICAL SOL WITH MATHEMATICA ===========================

St   = 1;
XiXi = sigma_a1^2;
yin  = 0;
vin  = 0;
u    = 1;

for j=1:length(P1.t)
    t = P1.t(j);
    % ---------------------------------------------------------------------
    St = 1;
    mean_Yp1(j) = t.*u+((-1)+exp(1).^((-1).*St.^(-1).*t)).*St.*(u+(-1).*vin)+(1/2).* ...
        St.^(-1).*((-2).*St.^2+exp(1).^((-1).*St.^(-1).*t).*(2.*St.^2+2.* ...
        St.*t+t.^2)).*(u+(-1).*vin).*XiXi+yin;
    mean_Vp1(j) = u+exp(1).^((-1).*St.^(-1).*t).*((-1).*u+vin)+(1/2).*exp(1).^((-1) ...
        .*St.^(-1).*t).*St.^(-2).*t.^2.*((-1).*u+vin).*XiXi;
    YpYp1(j) = (1/4).*(u+(-1).*vin).^2.*XiXi.*(4.*(St+(-1).*exp(1).^((-1).*St.^( ...
        -1).*t).*(St+t)).^2+3.*St.^(-2).*((-2).*St.^2+exp(1).^((-1).*St.^( ...
        -1).*t).*(2.*St.^2+2.*St.*t+t.^2)).^2.*XiXi);
    YpVp1(j) = (-1/4).*exp(1).^((-2).*St.^(-1).*t).*St.^(-3).*t.*(u+(-1).*vin) ...
        .^2.*XiXi.*(4.*St.^3+4.*St.^2.*t+6.*St.^2.*t.*XiXi+6.*St.*t.^2.* ...
        XiXi+3.*t.^3.*XiXi+(-2).*exp(1).^(St.^(-1).*t).*St.^2.*(2.*St+3.* ...
        t.*XiXi));
    VpVp1(j) = (1/4).*exp(1).^((-2).*St.^(-1).*t).*St.^(-4).*t.^2.*(u+(-1).*vin) ...
        .^2.*XiXi.*(4.*St.^2+3.*t.^2.*XiXi);
    % ---------------------------------------------------------------------
    St = 5;
    mean_Yp5(j) = t.*u+((-1)+exp(1).^((-1).*St.^(-1).*t)).*St.*(u+(-1).*vin)+(1/2).* ...
        St.^(-1).*((-2).*St.^2+exp(1).^((-1).*St.^(-1).*t).*(2.*St.^2+2.* ...
        St.*t+t.^2)).*(u+(-1).*vin).*XiXi+yin;
    mean_Vp5(j) = u+exp(1).^((-1).*St.^(-1).*t).*((-1).*u+vin)+(1/2).*exp(1).^((-1) ...
        .*St.^(-1).*t).*St.^(-2).*t.^2.*((-1).*u+vin).*XiXi;
    YpYp5(j) = (1/4).*(u+(-1).*vin).^2.*XiXi.*(4.*(St+(-1).*exp(1).^((-1).*St.^( ...
        -1).*t).*(St+t)).^2+3.*St.^(-2).*((-2).*St.^2+exp(1).^((-1).*St.^( ...
        -1).*t).*(2.*St.^2+2.*St.*t+t.^2)).^2.*XiXi);
    YpVp5(j) = (-1/4).*exp(1).^((-2).*St.^(-1).*t).*St.^(-3).*t.*(u+(-1).*vin) ...
        .^2.*XiXi.*(4.*St.^3+4.*St.^2.*t+6.*St.^2.*t.*XiXi+6.*St.*t.^2.* ...
        XiXi+3.*t.^3.*XiXi+(-2).*exp(1).^(St.^(-1).*t).*St.^2.*(2.*St+3.* ...
        t.*XiXi));
    VpVp5(j) = (1/4).*exp(1).^((-2).*St.^(-1).*t).*St.^(-4).*t.^2.*(u+(-1).*vin) ...
        .^2.*XiXi.*(4.*St.^2+3.*t.^2.*XiXi);
end

close all

subplot(221)
plot(P1.t,squeeze(P1.MoM.mean_xp(1,1,:))); hold on
plot(P1.t,mean_Yp1,'--')
plot(P1.t,squeeze(P1.MoM.mean_up(1,1,:))); hold on
plot(P1.t,mean_Vp1,'--')
xlabel('t')
ylabel('First moments')
title('St=1')

subplot(222)
plot(P1.t,squeeze(P1.MoM.xpxp(1,1,:))); hold on
plot(P1.t,YpYp1,'--')
plot(P1.t,squeeze(P1.MoM.xpup(1,1,:))); hold on
plot(P1.t,YpVp1,'--')
plot(P1.t,squeeze(P1.MoM.upup(1,1,:))); hold on
plot(P1.t,VpVp1,'--')
xlabel('t')
ylabel('Second moments')
title('St=1')

subplot(223)
plot(P1.t,squeeze(P5.MoM.mean_xp(1,1,:))); hold on
plot(P1.t,mean_Yp5,'--')
plot(P1.t,squeeze(P5.MoM.mean_up(1,1,:))); hold on
plot(P1.t,mean_Vp5,'--')
xlabel('t')
ylabel('First moments')
title('St=5')

subplot(224)
plot(P1.t,squeeze(P5.MoM.xpxp(1,1,:))); hold on
plot(P1.t,YpYp5,'--')
plot(P1.t,squeeze(P5.MoM.xpup(1,1,:))); hold on
plot(P1.t,YpVp5,'--')
plot(P1.t,squeeze(P5.MoM.upup(1,1,:))); hold on
plot(P1.t,VpVp5,'--')
xlabel('t')
ylabel('Second moments')
title('St=5')






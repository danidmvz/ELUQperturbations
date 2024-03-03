%% ========================================================================
close all
clear
clc
setup
label    = 'doubleGyre';
iFig     = 1;
saveFigs = 1;

% PARAMETERS --------------------------------------------------------------
timeMethod = 1;
flowType   = 3;
dataFlow.e = 0.25;
dataFlow.w = 2*pi/10;
dataFlow.A = 0.1;
f1Type     = 3; % Stokes random
dataf1     = [];
npx        = 201; % Numer of clouds along x
npy        = round(npx/2); % Numer of clouds along y
nt         = [251 2];
tLim       = [0 6];

% Forcing and initial condition
mean_a1  = 1;
sigma_a1 = 0.1;
cm3_a1   = 0;
lim1     = mean_a1-0.5*sqrt(12)*sigma_a1;
lim2     = mean_a1+0.5*sqrt(12)*sigma_a1;
cm4_a1   = (lim2-lim1)^4/80;
x        = linspace(0,2,npx)';
y        = linspace(0,1,npy)';
for i=1:npx
    for j=1:npy
        xp00(i,j) = x(i);
        yp00(i,j) = y(j);
        up00(i,j) = 0;
        vp00(i,j) = 0;
    end
end


%% MC-PSIC ================================================================

% -------------------------------------------------------------------------
taup = 1;
ns   = 100; % Numer of samples per cloud in MC-PSIC
% a1 = unifrnd(mean_a1-0.5*sqrt(12)*sigma_a1,mean_a1+0.5*sqrt(12)*sigma_a1,ns,1);
a1 = linspace(mean_a1-0.5*sqrt(12)*sigma_a1,mean_a1+0.5*sqrt(12)*sigma_a1,ns)';
for i=1:npx
    for j=1:npy
        for k=1:ns
            xp0(i,j,k)       = xp00(i,j);
            yp0(i,j,k)       = yp00(i,j);
            up0(i,j,k)       = up00(i,j);
            vp0(i,j,k)       = vp00(i,j);
            dataf1.a1(i,j,k) = a1(k);
        end
    end
end
M1 = FUNC_solverMCPSIC2D_inertial(xp0,yp0,up0,vp0,tLim,nt,flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

% -------------------------------------------------------------------------
taup = 5;
M5 = FUNC_solverMCPSIC2D_inertial(xp0,yp0,up0,vp0,tLim,nt,flowType,dataFlow,f1Type,dataf1,timeMethod,taup);


%% MoP ====================================================================

% -------------------------------------------------------------------------
taup = 1;
clear xp0 yp0 up0 vp0
for i=1:npx
    for j=1:npy
        xp0(i,j) = xp00(i,j);
        yp0(i,j) = yp00(i,j);
        up0(i,j) = up00(i,j);
        vp0(i,j) = vp00(i,j);
    end
end
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

%% SPARSER ================================================================

% -------------------------------------------------------------------------
taup     = 1;
mean_xp0 = xp00;
mean_yp0 = yp00;
mean_up0 = up00;
mean_vp0 = vp00;
mean_a   = mean_a1;
aa       = sigma_a1^2;
S11 = FUNC_solverSPARSER2D_inertial(mean_xp0,mean_yp0,mean_up0,mean_vp0, ...
    mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,1,taup);
% S12 = FUNC_solverSPARSER2D_inertial(mean_xp0,mean_yp0,mean_up0,mean_vp0, ...
%     mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,2,taup);
% S13 = FUNC_solverSPARSER2D_inertial(mean_xp0,mean_yp0,mean_up0,mean_vp0, ...
%     mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,3,taup);
% S14 = FUNC_solverSPARSER2D_inertial(mean_xp0,mean_yp0,mean_up0,mean_vp0, ...
%     mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,4,taup);
% S15 = FUNC_solverSPARSER2D_inertial(mean_xp0,mean_yp0,mean_up0,mean_vp0, ...
%     mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,5,taup);

% -------------------------------------------------------------------------
taup = 5;
S51 = FUNC_solverSPARSER2D_inertial(mean_xp0,mean_yp0,mean_up0,mean_vp0, ...
    mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,1,taup);
% S52 = FUNC_solverSPARSER2D_inertial(mean_xp0,mean_yp0,mean_up0,mean_vp0, ...
%     mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,2,taup);
% S53 = FUNC_solverSPARSER2D_inertial(mean_xp0,mean_yp0,mean_up0,mean_vp0, ...
%     mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,3,taup);
% S54 = FUNC_solverSPARSER2D_inertial(mean_xp0,mean_yp0,mean_up0,mean_vp0, ...
%     mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,4,taup);
% S55 = FUNC_solverSPARSER2D_inertial(mean_xp0,mean_yp0,mean_up0,mean_vp0, ...
%     mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,5,taup);


%% PLOTS: CONTOURS ========================================================

close all
levels = 500;

j  = nt(1);
xx1 = 1.38;
yy1 = 0.3;
xx2 = 0.9;
yy2 = 0.74;
[~,ix1] = min(abs(x-xx1));
[~,iy1] = min(abs(y-yy1));
[~,ix2] = min(abs(x-xx2));
[~,iy2] = min(abs(y-yy2));

% fig = figure('units','inches','outerposition',[1 1 8+1 5+1],'color','w'); iFig = fig.Number;
figure(iFig)
contourf(x,y,(M1.maxEigKxp(:,:,j)').^0.5,levels,'edgecolor','none'); hold on
scatter(xx1,yy1,50,'w','d','filled')
shading interp
colormap jet
cb = colorbar;
cb.FontSize = 0.8*fs;
cb.TickLabelInterpreter = 'latex';
cTH = get(cb,'Title');
% xlim([min(x) max(x)])
% ylim([min(y) max(y)])
caxis([0 0.1]);
cb.Ticks = 0:0.02:0.1;
axis equal
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
xlabel('$y_1$','interpreter','latex','fontsize',fs)
ylabel('$y_2$','interpreter','latex','fontsize',fs)
title('$\hat{y}_{t_0}^t$','interpreter','latex','fontsize',fs)
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'dG_lambda_y_St1';
    savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end

% fig = figure('units','inches','outerposition',[1 1 8+1 5+1],'color','w'); iFig = fig.Number;
figure(iFig)
contourf(x,y,(M5.maxEigKxp(:,:,j)').^0.5,levels,'edgecolor','none'); hold on
shading interp
colormap jet
cb = colorbar;
cb.FontSize = 0.8*fs;
cb.TickLabelInterpreter = 'latex';
cTH = get(cb,'Title');
% xlim([min(x) max(x)])
% ylim([min(y) max(y)])
caxis([0 0.1]);
cb.Ticks = 0:0.02:0.1;
axis equal
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
xlabel('$y_1$','interpreter','latex','fontsize',fs)
ylabel('$y_2$','interpreter','latex','fontsize',fs)
title('$\hat{y}_{t_0}^t$','interpreter','latex','fontsize',fs)
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'dG_lambda_y_St5';
    savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end


% fig = figure('units','inches','outerposition',[1 1 8+1 5+1],'color','w'); iFig = fig.Number;
figure(iFig)
contourf(x,y,(M1.maxEigKup(:,:,j)').^0.5,levels,'edgecolor','none'); hold on
shading interp
colormap jet
cb = colorbar;
cb.FontSize = 0.8*fs;
cb.TickLabelInterpreter = 'latex';
cTH = get(cb,'Title');
xlim([min(x) max(x)])
ylim([min(y) max(y)])
caxis([0 0.05]);
cb.Ticks = 0:0.01:0.05;
axis equal
hold off
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
xlabel('$y_1$','interpreter','latex','fontsize',fs)
ylabel('$y_2$','interpreter','latex','fontsize',fs)
title('$\hat{v}_{t_0}^t$','interpreter','latex','fontsize',fs)
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'dG_lambda_v_St1';
    savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end


% fig = figure('units','inches','outerposition',[1 1 8+1 5+1],'color','w'); iFig = fig.Number;
figure(iFig)
contourf(x,y,(M5.maxEigKup(:,:,j)').^0.5,levels,'edgecolor','none'); hold on
shading interp
colormap jet
cb = colorbar;
cb.FontSize = 0.8*fs;
cb.TickLabelInterpreter = 'latex';
cTH = get(cb,'Title');
xlim([min(x) max(x)])
ylim([min(y) max(y)])
caxis([0 0.05]);
cb.Ticks = 0:0.01:0.05;
axis equal
hold off
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
xlabel('$y_1$','interpreter','latex','fontsize',fs)
ylabel('$y_2$','interpreter','latex','fontsize',fs)
title('$\hat{v}_{t_0}^t$','interpreter','latex','fontsize',fs)
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'dG_lambda_v_St5';
    savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end


%% PLOTS: PATH ============================================================

clear xp0 yp0 up0 vp0 dataf1
taup = 1.5;
ns   = 100; % Numer of samples per cloud in MC-PSIC
a1 = linspace(mean_a1-0.5*sqrt(12)*sigma_a1,mean_a1+0.5*sqrt(12)*sigma_a1,ns)';
for i=1:npx
    for j=1:npy
        for k=1:ns
            xp0(1,1,k)       = xx1;
            yp0(1,1,k)       = yy1;
            up0(1,1,k)       = 0;
            vp0(1,1,k)       = 0;
            dataf1.a1(1,1,k) = a1(k);
        end
    end
end
m15 = FUNC_solverMCPSIC2D_inertial(xp0,yp0,up0,vp0,tLim,nt,flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

clear xp0 yp0 up0 vp0 dataf1
taup = 0.5;
ns   = 100; % Numer of samples per cloud in MC-PSIC
a1 = linspace(mean_a1-0.5*sqrt(12)*sigma_a1,mean_a1+0.5*sqrt(12)*sigma_a1,ns)';
for i=1:npx
    for j=1:npy
        for k=1:ns
            xp0(1,1,k)       = xx1;
            yp0(1,1,k)       = yy1;
            up0(1,1,k)       = 0;
            vp0(1,1,k)       = 0;
            dataf1.a1(1,1,k) = a1(k);
        end
    end
end
m05 = FUNC_solverMCPSIC2D_inertial(xp0,yp0,up0,vp0,tLim,nt,flowType,dataFlow,f1Type,dataf1,timeMethod,taup);


close all
fig = figure('units','inches','outerposition',[1 1 0.8*(15+1) 0.8*(6+1)],'color','w'); iFig = fig.Number;
p1(2) = plot(squeeze(m05.mean_xp(1,1,:))   ,squeeze(m05.mean_yp(1,1,:))   ,'-k','color',color14,'linewidth',0.5*lw); hold on
p1(1) = plot(squeeze(M1.mean_xp(ix1,iy1,:)),squeeze(M1.mean_yp(ix1,iy1,:)),'-k','color',color14,'linewidth',0.5*lw); hold on
p1(3) = plot(squeeze(m15.mean_xp(1,1,:))   ,squeeze(m15.mean_yp(1,1,:))   ,'-k','color',color14,'linewidth',0.5*lw); hold on
jj    = round(linspace(1,nt(1),15));
for j=jj
    
    [X,Y] = FUNC_ellipse(squeeze(P1.mean_xp(ix1,iy1,j)),squeeze(P1.mean_yp(ix1,iy1,j)), ...
        squeeze(P1.xpxp(ix1,iy1,j)),squeeze(P1.ypyp(ix1,iy1,j)),squeeze(P1.xpyp(ix1,iy1,j)),1e2);
    fill(X,Y,color9,'linestyle','none','facealpha',0.2);
    plot(X,Y,'-','color',color9,'linewidth',lw)
    plot(squeeze(M1.xp(ix1,iy1,:,j)),squeeze(M1.yp(ix1,iy1,:,j)),'r','color',color10,'linewidth',0.5*lw)
    plot(squeeze(M1.mean_xp(ix1,iy1,j)),squeeze(M1.mean_yp(ix1,iy1,j)),'.k','linewidth',lw); hold on
    
    [X,Y] = FUNC_ellipse(squeeze(m05.mean_xp(1,1,j)),squeeze(m05.mean_yp(1,1,j)), ...
        squeeze(m05.xpxp(1,1,j)),squeeze(m05.ypyp(1,1,j)),squeeze(m05.xpyp(1,1,j)),1e2);
    fill(X,Y,color9,'linestyle','none','facealpha',0.2);
    plot(X,Y,'-','color',color9,'linewidth',lw)
    plot(squeeze(m05.xp(1,1,:,j)),squeeze(m05.yp(1,1,:,j)),'r','color',color10,'linewidth',0.5*lw)
    plot(squeeze(m05.mean_xp(1,1,j)),squeeze(m05.mean_yp(1,1,j)),'.k','linewidth',lw); hold on
    
    [X,Y] = FUNC_ellipse(squeeze(m15.mean_xp(1,1,j)),squeeze(m15.mean_yp(1,1,j)), ...
        squeeze(m15.xpxp(1,1,j)),squeeze(m15.ypyp(1,1,j)),squeeze(m15.xpyp(1,1,j)),1e2);
    fill(X,Y,color9,'linestyle','none','facealpha',0.2);
    plot(X,Y,'-','color',color9,'linewidth',lw)
    plot(squeeze(m15.xp(1,1,:,j)),squeeze(m15.yp(1,1,:,j)),'r','color',color10,'linewidth',0.5*lw)
    plot(squeeze(m15.mean_xp(1,1,j)),squeeze(m15.mean_yp(1,1,j)),'.k','linewidth',lw); hold on
    
    drawnow
end
scatter(xx1,yy1,50,'k','d','filled')
text(2.13+0.03,0.2,'$St=1.5$','interpreter','latex','fontsize',fs)
text(2.1+0.03,0.5,'$St=1$','interpreter','latex','fontsize',fs)
text(2.08+0.03,1,'$St=0.5$','interpreter','latex','fontsize',fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
xlabel('$y_1$','interpreter','latex','fontsize',fs)
ylabel('$y_2$','interpreter','latex','fontsize',fs)
xlim([1.3 2.3])
% ax.YAxis.TickValues = 0:0.1:0.7;
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'dG_ellipses';
    savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg


%% PLOTS: MOMENTS =========================================================

nn = round(0.5*nt(1));
i = ix1;
j = iy1;

close all

figure(iFig) % close all
p1(1) = plot(M1.t,squeeze(M1.xpxp(i,j,:).^0.5),'s-','linewidth',lw,'color',color14,'MarkerIndices',round(0.1*nn):nn:nt(1),'MarkerFaceColor',color14); hold on
p1(2) = plot(P1.t,squeeze(P1.xpxp(i,j,:).^0.5),'>--','linewidth',lw,'color',color8,'MarkerIndices',round(0.3*nn):nn:nt(1),'MarkerFaceColor',color8); hold on
p1(3) = plot(P1.t,squeeze(P1.MoM.xpxp(i,j,:).^0.5),'o-.','linewidth',lw,'color',color10,'MarkerIndices',round(0.6*nn):nn:nt(1),'MarkerFaceColor',color10); hold on
p1(4) = plot(S11.t,squeeze(S11.xpxp(i,j,:).^0.5),'hexagram:','linewidth',lw,'color',color9,'MarkerIndices',round(0.8*nn):nn:nt(1),'MarkerFaceColor',color9); hold on
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylim([0 0.04])
yyaxis right
p1(5) = plot(M1.t,squeeze(M1.ypyp(i,j,:).^0.5),'d-','linewidth',lw,'color',color15,'MarkerIndices',round(0.1*nn):nn:nt(1),'MarkerFaceColor',color15); hold on
p1(6) = plot(P1.t,squeeze(P1.ypyp(i,j,:).^0.5),'^--','linewidth',lw,'color',color1,'MarkerIndices',round(0.3*nn):nn:nt(1),'MarkerFaceColor',color1); hold on
p1(7) = plot(P1.t,squeeze(P1.MoM.ypyp(i,j,:).^0.5),'*-.','linewidth',lw,'color',color2,'MarkerIndices',round(0.6*nn):nn:nt(1),'MarkerFaceColor',color2); hold on
p1(8) = plot(S11.t,squeeze(S11.ypyp(i,j,:).^0.5),'pentagram:','linewidth',lw,'color',color3,'MarkerIndices',round(0.8*nn):nn:nt(1),'MarkerFaceColor',color3); hold on
ylim([0 0.14])
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ax.YTick = 0:0.02:0.2;
ylabel('$\sigma_{Y_2}$','interpreter','latex','fontsize',fs)
xlabel('$t$','interpreter','latex','fontsize',fs)
yyaxis left
ylabel('$\sigma_{Y_1}$','interpreter','latex','fontsize',fs)
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
leg1{1} = '$\sigma_{Y_1}$, MC-PSIC';
leg1{2} = '$\sigma_{Y_1}$, MC-MoP';
leg1{3} = '$\sigma_{Y_1}$, MoM-MoP';
leg1{4} = '$\sigma_{Y_1}$, SPARSE-R';
leg1{5} = '$\sigma_{Y_2}$, MC-PSIC';
leg1{6} = '$\sigma_{Y_2}$, MC-MoP';
leg1{7} = '$\sigma_{Y_2}$, MoM-MoP';
leg1{8} = '$\sigma_{Y_2}$, SPARSE-R';
leg = legend(p1,leg1,'interpreter','latex','fontsize',0.68*fs,'location','north');
leg.Box = 'off';
leg.NumColumns = 2;
ax.XTick = 0:1:10;
ax.YTick = 0:0.01:0.2;
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'dG_cm2_xx';
    savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg


figure(iFig) % close all
p1(1) = plot(M1.t,squeeze(M1.upup(i,j,:).^0.5),'s-','linewidth',lw,'color',color14,'MarkerIndices',round(0.1*nn):nn:nt(1),'MarkerFaceColor',color14); hold on
p1(2) = plot(P1.t,squeeze(P1.upup(i,j,:).^0.5),'>--','linewidth',lw,'color',color8,'MarkerIndices',round(0.3*nn):nn:nt(1),'MarkerFaceColor',color8); hold on
p1(3) = plot(P1.t,squeeze(P1.MoM.upup(i,j,:).^0.5),'o-.','linewidth',lw,'color',color10,'MarkerIndices',round(0.6*nn):nn:nt(1),'MarkerFaceColor',color10); hold on
p1(4) = plot(S11.t,squeeze(S11.upup(i,j,:).^0.5),'hexagram:','linewidth',lw,'color',color9,'MarkerIndices',round(0.8*nn):nn:nt(1),'MarkerFaceColor',color9); hold on
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
ylim([0 0.04])
yyaxis right
p1(5) = plot(M1.t,squeeze(M1.vpvp(i,j,:).^0.5),'d-','linewidth',lw,'color',color15,'MarkerIndices',round(0.1*nn):nn:nt(1),'MarkerFaceColor',color15); hold on
p1(6) = plot(P1.t,squeeze(P1.vpvp(i,j,:).^0.5),'^--','linewidth',lw,'color',color1,'MarkerIndices',round(0.3*nn):nn:nt(1),'MarkerFaceColor',color1); hold on
p1(7) = plot(P1.t,squeeze(P1.MoM.vpvp(i,j,:).^0.5),'*-.','linewidth',lw,'color',color2,'MarkerIndices',round(0.6*nn):nn:nt(1),'MarkerFaceColor',color2); hold on
p1(8) = plot(S11.t,squeeze(S11.vpvp(i,j,:).^0.5),'pentagram:','linewidth',lw,'color',color3,'MarkerIndices',round(0.8*nn):nn:nt(1),'MarkerFaceColor',color3); hold on
ylim([0 0.08])
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
% ax.YTick = 0:0.02:0.2;
ylabel('$\sigma_{V_2}$','interpreter','latex','fontsize',fs)
xlabel('$t$','interpreter','latex','fontsize',fs)
yyaxis left
ylabel('$\sigma_{V_1}$','interpreter','latex','fontsize',fs)
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
leg1{1} = '$\sigma_{V_1}$, MC-PSIC';
leg1{2} = '$\sigma_{V_1}$, MC-MoP';
leg1{3} = '$\sigma_{V_1}$, MoM-MoP';
leg1{4} = '$\sigma_{V_1}$, SPARSE-R';
leg1{5} = '$\sigma_{V_2}$, MC-PSIC';
leg1{6} = '$\sigma_{V_2}$, MC-MoP';
leg1{7} = '$\sigma_{V_2}$, MoM-MoP';
leg1{8} = '$\sigma_{V_2}$, SPARSE-R';
leg = legend(p1,leg1,'interpreter','latex','fontsize',0.68*fs,'location','north');
leg.Box = 'off';
leg.NumColumns = 2;
ax.XTick = 0:1:10;
% ax.YTick = 0:0.01:0.2;
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'dG_cm2_uu';
    savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg


%% PLOTS: ORDERS ==========================================================

i = ix1;
j = iy1;

% -------------------------------------------------------------------------
figure(iFig) % close all
p1(1) = plot(P1.t,squeeze(P1.o0.xp(i,j,:)),'linewidth',lw,'color',color14); hold on
p1(2) = plot(P1.t,squeeze(P1.o1.xp(i,j,:)),'--','linewidth',lw,'color',color14); hold on
p1(3) = plot(P1.t,squeeze(P1.o2.xp(i,j,:)),'-.','linewidth',lw,'color',color14); hold on
p1(1) = plot(P1.t,squeeze(P1.o0.yp(i,j,:)),'linewidth',lw,'color',color8); hold on
p1(2) = plot(P1.t,squeeze(P1.o1.yp(i,j,:)),'--','linewidth',lw,'color',color8); hold on
p1(3) = plot(P1.t,squeeze(P1.o2.yp(i,j,:)),'-.','linewidth',lw,'color',color8); hold on

p1(1) = plot(P1.t,squeeze(P1.o0.up(i,j,:)),'linewidth',lw,'color',color9); hold on
p1(2) = plot(P1.t,squeeze(P1.o1.up(i,j,:)),'--','linewidth',lw,'color',color9); hold on
p1(3) = plot(P1.t,squeeze(P1.o2.up(i,j,:)),'-.','linewidth',lw,'color',color9); hold on
p1(1) = plot(P1.t,squeeze(P1.o0.vp(i,j,:)),'linewidth',lw,'color',color10); hold on
p1(2) = plot(P1.t,squeeze(P1.o1.vp(i,j,:)),'--','linewidth',lw,'color',color10); hold on
p1(3) = plot(P1.t,squeeze(P1.o2.vp(i,j,:)),'-.','linewidth',lw,'color',color10); hold on

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








%% ========================================================================
close all
clear
clc
setup
label    = 'vKVS';
iFig     = 1;
saveFigs = 1;

% PARAMETERS --------------------------------------------------------------
timeMethod = 1;
flowType   = 2;
dataFlow.Uref = 1;
f1Type     = 3; % Stokes random
dataf1     = [];
npx        = 201; % Numer of clouds along x
npy        = round(npx/2); % Numer of clouds along y
nt         = [151 2];
tLim       = [0 2.5];

% Forcing and initial condition
sigma_a1 = 0.1;
cm3_a1   = 0;
cm4_a1   = 0;
mean_a1  = 1;
x        = linspace(-1.5,5  ,npx)';
y        = linspace(-1.5,1.5,npy)';
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
% taup = 5;
% M5 = FUNC_solverMCPSIC2D_inertial(xp0,yp0,up0,vp0,tLim,nt,flowType,dataFlow,f1Type,dataf1,timeMethod,taup);


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
dataf1.a1       = a1;
dataf1.a1a1     = sigma_a1^2;
dataf1.a1a1a1   = cm3_a1;
dataf1.a1a1a1a1 = cm4_a1;
P1 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,tLim,nt, ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup);

% -------------------------------------------------------------------------
% taup = 5;
% P5 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,tLim,nt, ...
%     flowType,dataFlow,f1Type,dataf1,timeMethod,taup);


%% ANIMATION ==============================================================

levels = 100;
% fs     = 20;

close all
% cd('videos')
% name = 'dG_videoComparison_timeLegend';
% wobj = VideoWriter([name '.mp4'],'MPEG-4'); % close(wobj);
% wobj.FrameRate = 24; % frames per second (video speed)
% open(wobj)
% cont=1;
figure('units','inches','outerposition',[1 1 8+1 8+1],'color','w')
for j=1:10:nt(1)
    
    subplot(221)
    contourf(x,y,(P1.MoM.maxEigKxp(:,:,j)').^0.5,levels,'edgecolor','none'); hold on
    shading interp
    colormap jet
    cb = colorbar;
    cb.FontSize = 0.6*fs;
    cb.TickLabelInterpreter = 'latex';
    cTH = get(cb,'Title');
    xlim([min(x) max(x)])
    ylim([min(y) max(y)])
    %     caxis([min(min(P1.MoM.maxEigKxp(:,:,j).^0.5)) max(max(P1.MoM.maxEigKxp(:,:,j).^0.5))]);
    axis equal
    hold off
    xlabel('$x$','interpreter','latex','fontsize',fs)
    ylabel('$y$','interpreter','latex','fontsize',fs)
    title('MoP, $\sqrt{\textrm{max}(\lambda_y)}$','interpreter','latex','fontsize',fs)
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize=0.6*fs;
    
    subplot(222)
    contourf(x,y,(P1.MoM.maxEigKup(:,:,j)').^0.5,levels,'edgecolor','none'); hold on
    shading interp
    colormap jet
    cb = colorbar;
    cb.FontSize = 0.6*fs;
    cb.TickLabelInterpreter = 'latex';
    cTH = get(cb,'Title');
    xlim([min(x) max(x)])
    ylim([min(y) max(y)])
    %     caxis([min(min(P1.MoM.maxEigKup(:,:,j).^0.5)) max(max(P1.MoM.maxEigKup(:,:,j).^0.5))]);
    axis equal
    hold off
    xlabel('$x$','interpreter','latex','fontsize',fs)
    ylabel('$y$','interpreter','latex','fontsize',fs)
    title('MoP, $\sqrt{\textrm{max}(\lambda_v)}$','interpreter','latex','fontsize',fs)
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize=0.6*fs;
    
    
    subplot(223)
    contourf(x,y,(M1.maxEigKxp(:,:,j)').^0.5,levels,'edgecolor','none'); hold on
    shading interp
    colormap jet
    cb = colorbar;
    cb.FontSize = 0.6*fs;
    cb.TickLabelInterpreter = 'latex';
    cTH = get(cb,'Title');
    xlim([min(x) max(x)])
    ylim([min(y) max(y)])
    %     caxis([min(min(P1.MoM.maxEigKxp(:,:,j).^0.5)) max(max(P1.MoM.maxEigKxp(:,:,j).^0.5))]);
    axis equal
    hold off
    xlabel('$x$','interpreter','latex','fontsize',fs)
    ylabel('$y$','interpreter','latex','fontsize',fs)
    title('MC-PSIC, $\sqrt{\textrm{max}(\lambda_y)}$','interpreter','latex','fontsize',fs)
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize=0.6*fs;
    
    subplot(224)
    contourf(x,y,(M1.maxEigKup(:,:,j)').^0.5,levels,'edgecolor','none'); hold on
    shading interp
    colormap jet
    cb = colorbar;
    cb.FontSize = 0.6*fs;
    cb.TickLabelInterpreter = 'latex';
    cTH = get(cb,'Title');
    xlim([min(x) max(x)])
    ylim([min(y) max(y)])
    %     caxis([min(min(P1.MoM.maxEigKup(:,:,j).^0.5)) max(max(P1.MoM.maxEigKup(:,:,j).^0.5))]);
    axis equal
    hold off
    xlabel('$x$','interpreter','latex','fontsize',fs)
    ylabel('$y$','interpreter','latex','fontsize',fs)
    title('MC-PSIC, $\sqrt{\textrm{max}(\lambda_v)}$','interpreter','latex','fontsize',fs)
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize=0.6*fs;
    
    drawnow
    
    j
end



%% PLOTS: CONTOURS ========================================================

close all
levels = 200;

j  = nt(1);
xx1 = 1.38;
yy1 = 0.3;
xx2 = 0.9;
yy2 = 0.74;
[~,ix1] = min(abs(x-xx1));
[~,iy1] = min(abs(y-yy1));
[~,ix2] = min(abs(x-xx2));
[~,iy2] = min(abs(y-yy2));

fig = figure('units','inches','outerposition',[1 1 8+1 5+1],'color','w'); iFig = fig.Number;
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
% caxis([0 0.1]);
% cb.Ticks = 0:0.02:0.1;
axis equal
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
xlabel('$x$','interpreter','latex','fontsize',fs)
ylabel('$y$','interpreter','latex','fontsize',fs)
title('$\hat{y}_{t_0}^t$','interpreter','latex','fontsize',0.8*fs)
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'vKVSlambda_y_St1';
    %     savefig(iFig,thisName)
    %     print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end

% fig = figure('units','inches','outerposition',[1 1 8+1 5+1],'color','w'); iFig = fig.Number;
% contourf(x,y,(M5.maxEigKxp(:,:,j)').^0.5,levels,'edgecolor','none'); hold on
% shading interp
% colormap jet
% cb = colorbar;
% cb.FontSize = 0.8*fs;
% cb.TickLabelInterpreter = 'latex';
% cTH = get(cb,'Title');
% % xlim([min(x) max(x)])
% % ylim([min(y) max(y)])
% caxis([0 0.1]);
% cb.Ticks = 0:0.02:0.1;
% axis equal
% ax = gca;
% ax.TickLabelInterpreter = 'latex';
% ax.FontSize=0.8*fs;
% xlabel('$x$','interpreter','latex','fontsize',fs)
% ylabel('$y$','interpreter','latex','fontsize',fs)
% title('$\hat{y}_{t_0}^t$','interpreter','latex','fontsize',0.8*fs)
% if saveFigs == 1
%     cd([figuresDir '_' label]);
%     thisName = 'vKVSlambda_y_St5';
%     %     savefig(iFig,thisName)
%     %     print(iFig,thisName,'-depsc','-painters')
%     print(iFig,thisName,'-dpng','-r800','-painters')
%     cd ..;
% end


fig = figure('units','inches','outerposition',[1 1 8+1 5+1],'color','w'); iFig = fig.Number;
contourf(x,y,(M1.maxEigKup(:,:,j)').^0.5,levels,'edgecolor','none'); hold on
shading interp
colormap jet
cb = colorbar;
cb.FontSize = 0.8*fs;
cb.TickLabelInterpreter = 'latex';
cTH = get(cb,'Title');
xlim([min(x) max(x)])
ylim([min(y) max(y)])
% caxis([0 0.05]);
% cb.Ticks = 0:0.01:0.05;
axis equal
hold off
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
xlabel('$x$','interpreter','latex','fontsize',fs)
ylabel('$y$','interpreter','latex','fontsize',fs)
title('$\hat{v}_{t_0}^t$','interpreter','latex','fontsize',0.8*fs)
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'vKVSlambda_v_St1';
    %     savefig(iFig,thisName)
    %     print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end


% fig = figure('units','inches','outerposition',[1 1 8+1 5+1],'color','w'); iFig = fig.Number;
% contourf(x,y,(M5.maxEigKup(:,:,j)').^0.5,levels,'edgecolor','none'); hold on
% shading interp
% colormap jet
% cb = colorbar;
% cb.FontSize = 0.8*fs;
% cb.TickLabelInterpreter = 'latex';
% cTH = get(cb,'Title');
% xlim([min(x) max(x)])
% ylim([min(y) max(y)])
% caxis([0 0.05]);
% cb.Ticks = 0:0.01:0.05;
% axis equal
% hold off
% ax = gca;
% ax.TickLabelInterpreter = 'latex';
% ax.FontSize=0.8*fs;
% xlabel('$x$','interpreter','latex','fontsize',fs)
% ylabel('$y$','interpreter','latex','fontsize',fs)
% title('$\hat{v}_{t_0}^t$','interpreter','latex','fontsize',0.8*fs)
% if saveFigs == 1
%     cd([figuresDir '_' label]);
%     thisName = 'vKVSlambda_v_St5';
%     %     savefig(iFig,thisName)
%     %     print(iFig,thisName,'-depsc','-painters')
%     print(iFig,thisName,'-dpng','-r800','-painters')
%     cd ..;
% end


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
fig = figure('units','inches','outerposition',[1 1 15+1 6+1],'color','w'); iFig = fig.Number;
p1(1) = plot(squeeze(M1.mean_xp(ix1,iy1,:)),squeeze(M1.mean_yp(ix1,iy1,:)),'-k','linewidth',0.5*lw); hold on
p1(2) = plot(squeeze(m05.mean_xp(1,1,:)),squeeze(m05.mean_yp(1,1,:)),'-k','color',color14,'linewidth',0.5*lw); hold on
p1(3) = plot(squeeze(m15.mean_xp(1,1,:)),squeeze(m15.mean_yp(1,1,:)),'-k','color',color14,'linewidth',0.5*lw); hold on
for j=1:10:nt(1)
    
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
text(2.13,0.2,'$St=1.5$','interpreter','latex','fontsize',0.8*fs)
text(2.1,0.5,'$St=1$','interpreter','latex','fontsize',0.8*fs)
text(2.08,1,'$St=0.5$','interpreter','latex','fontsize',0.8*fs)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
xlabel('$x$','interpreter','latex','fontsize',fs)
ylabel('$y$','interpreter','latex','fontsize',fs)
% ax.YAxis.TickValues = 0:0.1:0.7;
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'vKVSellipses';
    %     savefig(iFig,thisName)
    %     print(iFig,thisName,'-depsc','-painters')
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
p1(2) = plot(P1.t,squeeze(P1.xpxp(i,j,:).^0.5),'>--','linewidth',lw,'color',color8,'MarkerIndices',round(0.4*nn):nn:nt(1),'MarkerFaceColor',color8); hold on
p1(3) = plot(P1.t,squeeze(P1.MoM.xpxp(i,j,:).^0.5),'o-.','linewidth',lw,'color',color10,'MarkerIndices',round(0.7*nn):nn:nt(1),'MarkerFaceColor',color10); hold on
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
% ylim([0 0.04])
yyaxis right
p1(4) = plot(M1.t,squeeze(M1.ypyp(i,j,:).^0.5),'d-','linewidth',lw,'color',color15,'MarkerIndices',round(0.1*nn):nn:nt(1),'MarkerFaceColor',color15); hold on
p1(5) = plot(P1.t,squeeze(P1.ypyp(i,j,:).^0.5),'^--','linewidth',lw,'color',color3,'MarkerIndices',round(0.4*nn):nn:nt(1),'MarkerFaceColor',color3); hold on
p1(6) = plot(P1.t,squeeze(P1.MoM.ypyp(i,j,:).^0.5),'*-.','linewidth',lw,'color',color9,'MarkerIndices',round(0.7*nn):nn:nt(1),'MarkerFaceColor',color9); hold on
% ylim([0 0.14])
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
% ax.YTick = 0:0.02:0.2;
ylabel('$\sigma_{Y_2}$','interpreter','latex','fontsize',fs)
xlabel('$t$','interpreter','latex','fontsize',fs)
yyaxis left
ylabel('$\sigma_{Y_1}$','interpreter','latex','fontsize',fs)
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
leg1{1} = '$\sigma_{Y_1}$, MC-PSIC';
leg1{2} = '$\sigma_{Y_1}$, MC-MoP';
leg1{3} = '$\sigma_{Y_1}$, MoM-MoP';
leg1{4} = '$\sigma_{Y_2}$, MC-PSIC';
leg1{5} = '$\sigma_{Y_2}$, MC-MoP';
leg1{6} = '$\sigma_{Y_2}$, MoM-MoP';
leg = legend(p1,leg1,'interpreter','latex','fontsize',0.68*fs,'location','north');
leg.Box = 'off';
leg.NumColumns = 2;
% ax.XTick = 0:1:10;
% ax.YTick = 0:0.01:0.2;
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'vKVScm2_xx';
    %     savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg


figure(iFig) % close all
p1(1) = plot(M1.t,squeeze(M1.upup(i,j,:).^0.5),'s-','linewidth',lw,'color',color14,'MarkerIndices',round(0.1*nn):nn:nt(1),'MarkerFaceColor',color14); hold on
p1(2) = plot(P1.t,squeeze(P1.upup(i,j,:).^0.5),'>--','linewidth',lw,'color',color8,'MarkerIndices',round(0.4*nn):nn:nt(1),'MarkerFaceColor',color8); hold on
p1(3) = plot(P1.t,squeeze(P1.MoM.upup(i,j,:).^0.5),'o-.','linewidth',lw,'color',color10,'MarkerIndices',round(0.7*nn):nn:nt(1),'MarkerFaceColor',color10); hold on
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;
% ylim([0 0.03])
yyaxis right
p1(4) = plot(M1.t,squeeze(M1.vpvp(i,j,:).^0.5),'d-','linewidth',lw,'color',color15,'MarkerIndices',round(0.1*nn):nn:nt(1),'MarkerFaceColor',color15); hold on
p1(5) = plot(P1.t,squeeze(P1.vpvp(i,j,:).^0.5),'^--','linewidth',lw,'color',color3,'MarkerIndices',round(0.4*nn):nn:nt(1),'MarkerFaceColor',color3); hold on
p1(6) = plot(P1.t,squeeze(P1.MoM.vpvp(i,j,:).^0.5),'*-.','linewidth',lw,'color',color9,'MarkerIndices',round(0.7*nn):nn:nt(1),'MarkerFaceColor',color9); hold on
% ylim([0 0.06])
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
leg1{4} = '$\sigma_{V_2}$, MC-PSIC';
leg1{5} = '$\sigma_{V_2}$, MC-MoP';
leg1{6} = '$\sigma_{V_2}$, MoM-MoP';
leg = legend(p1,leg1,'interpreter','latex','fontsize',0.68*fs,'location','northeast');
leg.Box = 'off';
leg.NumColumns = 2;
% ax.XTick = 0:1:10;
% ax.YTick = 0:0.01:0.2;
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'vKVScm2_uu';
    %     savefig(iFig,thisName)
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
% ylim([-2 10])
yyaxis right
p1(4) = plot(P1.t,squeeze(P1.o0.up(i,j,:)),'linewidth',lw,'color',color8); hold on
p1(5) = plot(P1.t,squeeze(P1.o1.up(i,j,:)),'--','linewidth',lw,'color',color8); hold on
p1(6) = plot(P1.t,squeeze(P1.o2.up(i,j,:)),'-.','linewidth',lw,'color',color8); hold on
% ylim([-0.5 1])
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
% ax.XTick = 0:2:10;
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'UF_orders_St1';
    %     savefig(iFig,thisName)
    print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end
iFig = iFig+1; clear p1 leg1 leg








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
nt         = [151 5];
tLim       = [0 6];
taup       = 1;

% Forcing and initial condition
sigma_a1 = 0.1;
cm3_a1   = 0;
cm4_a1   = 0;
mean_a1  = 1;
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

ns = 100; % Numer of samples per cloud in MC-PSIC
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


%% MoP ====================================================================

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


%% PLOTS ==================================================================

i = round(0.75*npx);
j = round(0.75*npy);

close all
subplot(221)
plot(M1.t,squeeze(M1.xpxp(i,j,:).^0.5)); hold on
plot(P1.t,squeeze(P1.xpxp(i,j,:).^0.5),'--'); hold on
plot(P1.t,squeeze(P1.MoM.xpxp(i,j,:).^0.5),'-.'); hold on

subplot(222)
plot(M1.t,squeeze(M1.ypyp(i,j,:).^0.5)); hold on
plot(P1.t,squeeze(P1.ypyp(i,j,:).^0.5),'--'); hold on
plot(P1.t,squeeze(P1.MoM.ypyp(i,j,:).^0.5),'-.'); hold on

subplot(223)
plot(M1.t,squeeze(M1.upup(i,j,:).^0.5)); hold on
plot(P1.t,squeeze(P1.upup(i,j,:).^0.5),'--'); hold on
plot(P1.t,squeeze(P1.MoM.upup(i,j,:).^0.5),'-.'); hold on

subplot(224)
plot(M1.t,squeeze(M1.vpvp(i,j,:).^0.5)); hold on
plot(P1.t,squeeze(P1.vpvp(i,j,:).^0.5),'--'); hold on
plot(P1.t,squeeze(P1.MoM.vpvp(i,j,:).^0.5),'-.'); hold on



% fig = figure('units','inches','outerposition',[1 1 8+1 6+1],'color','w'); iFig = fig.Number;
% xx = M0s1.xp(:,:,:,j);
% yy = M0s1.yp(:,:,:,j);
% scatter(xx(:),yy(:),0.001,'k')
% axis equal
% hold off
% xlim([0 2])
% ylim([0 1])
% ax=gca;
% ax.TickLabelInterpreter = 'latex';
% ax.FontSize=0.6*fs;
% ax.XTick = [];
% ax.YTick = [];
% title('$St = 0.5, \ \sigma_{f}=0$','interpreter','latex','fontsize',fs)
% if saveFigs == 1
%     cd([figuresDir '_' label]);
%     thisName = 'dG_intro_part_St1_a0_T12w6';
%     %     savefig(iFig,thisName)
%     %     print(iFig,thisName,'-depsc','-painters')
%     print(iFig,thisName,'-dpng','-r600','-painters')
%     cd ..;
% end
% iFig = iFig+1; clear p1 leg1 leg


%% MAX EIG ================================================================

levels = 100;
% fs     = 20;

% close all
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
    caxis([min(min(P1.MoM.maxEigKxp(:,:,j).^0.5)) max(max(P1.MoM.maxEigKxp(:,:,j).^0.5))]);
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
    caxis([min(min(P1.MoM.maxEigKup(:,:,j).^0.5)) max(max(P1.MoM.maxEigKup(:,:,j).^0.5))]);
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
    caxis([min(min(P1.MoM.maxEigKxp(:,:,j).^0.5)) max(max(P1.MoM.maxEigKxp(:,:,j).^0.5))]);
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
    caxis([min(min(P1.MoM.maxEigKup(:,:,j).^0.5)) max(max(P1.MoM.maxEigKup(:,:,j).^0.5))]);
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

%% PATH OF PARTICLES WITH UQ ==============================================

setup
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

close all

fig = figure('units','inches','outerposition',[1 1 8+1 5+1],'color','w'); iFig = fig.Number;
contourf(x,y,(M1.maxEigKxp(:,:,j)').^0.5,levels,'edgecolor','none'); hold on
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
xlabel('$x$','interpreter','latex','fontsize',fs)
ylabel('$y$','interpreter','latex','fontsize',fs)
title('$\hat{y}_{t_0}^t$','interpreter','latex','fontsize',0.8*fs)
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'dG_lambda_y_St1';
    %     savefig(iFig,thisName)
    %     print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end


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
caxis([0 0.05]);
cb.Ticks = 0:0.01:0.05;
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
    thisName = 'dG_lambda_v_St1';
    %     savefig(iFig,thisName)
    %     print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end

    
fig = figure('units','inches','outerposition',[1 1 15+1 4+1],'color','w'); iFig = fig.Number;
plot(squeeze(M1.mean_xp(ix1,iy1,:)),squeeze(M1.mean_yp(ix1,iy1,:)),'k','linewidth',0.5*lw); hold on
plot(squeeze(M1.mean_xp(ix1,iy1,1)),squeeze(M1.mean_yp(ix1,iy1,1)),'.k','linewidth',lw,'MarkerFaceColor',color14); hold on
for j=1:10:nt(1)
    [X,Y] = FUNC_ellipse(squeeze(P1.mean_xp(ix1,iy1,j)),squeeze(P1.mean_yp(ix1,iy1,j)), ...
        squeeze(P1.xpxp(ix1,iy1,j)),squeeze(P1.ypyp(ix1,iy1,j)),squeeze(P1.xpyp(ix1,iy1,j)),1e2);
    fill(X,Y,color9,'linestyle','none','facealpha',0.2);
    plot(X,Y,'-','color',color9,'linewidth',lw)
    plot(squeeze(P1.xp(ix1,iy1,:,j)),squeeze(P1.yp(ix1,iy1,:,j)),'r','linewidth',0.5*lw)
    plot(squeeze(M1.mean_xp(ix1,iy1,j)),squeeze(M1.mean_yp(ix1,iy1,j)),'.k','linewidth',lw); hold on
    drawnow
end
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize=0.8*fs;    
xlabel('$x$','interpreter','latex','fontsize',fs)
ylabel('$y$','interpreter','latex','fontsize',fs)
ax.YAxis.TickValues = 0:0.1:0.7;
if saveFigs == 1
    cd([figuresDir '_' label]);
    thisName = 'dG_ellipse_St1';
    %     savefig(iFig,thisName)
    %     print(iFig,thisName,'-depsc','-painters')
    print(iFig,thisName,'-dpng','-r800','-painters')
    cd ..;
end


%% MOMENTS ================================================================

nn = 41;
i = ix1;
j = iy1;

close all
subplot(221)
plot(M1.t,squeeze(M1.xpxp(i,j,:).^0.5),'d-','linewidth',lw,'color',color15,'MarkerIndices',round(0.1*nn):nn:nt(1),'MarkerFaceColor',color15);
plot(P1.t,squeeze(P1.xpxp(i,j,:).^0.5),'--'); hold on
plot(P1.t,squeeze(P1.MoM.xpxp(i,j,:).^0.5),'-.'); hold on

subplot(222)
plot(M1.t,squeeze(M1.ypyp(i,j,:).^0.5)); hold on
plot(P1.t,squeeze(P1.ypyp(i,j,:).^0.5),'--'); hold on
plot(P1.t,squeeze(P1.MoM.ypyp(i,j,:).^0.5),'-.'); hold on

subplot(223)
plot(M1.t,squeeze(M1.upup(i,j,:).^0.5)); hold on
plot(P1.t,squeeze(P1.upup(i,j,:).^0.5),'--'); hold on
plot(P1.t,squeeze(P1.MoM.upup(i,j,:).^0.5),'-.'); hold on

subplot(224)
plot(M1.t,squeeze(M1.vpvp(i,j,:).^0.5)); hold on
plot(P1.t,squeeze(P1.vpvp(i,j,:).^0.5),'--'); hold on
plot(P1.t,squeeze(P1.MoM.vpvp(i,j,:).^0.5),'-.'); hold on







%% 

levels = 100;
fs     = 20;

close all
% cd('videos')
% name = 'dG_videoComparison_timeLegend';
% wobj = VideoWriter([name '.mp4'],'MPEG-4'); % close(wobj);
% wobj.FrameRate = 24; % frames per second (video speed)
% open(wobj)
% cont=1;
figure('units','inches','outerposition',[1 1 8+1 8+1],'color','w')
for j=1:10:nt(1)
    
    subplot(311)
    contourf(x,y,sqrt(P1.o0.up(:,:,j)'.^2+P1.o0.vp(:,:,j)'.^2),levels,'edgecolor','none'); hold on
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
    xlabel('$x_1$','interpreter','latex','fontsize',fs)
    ylabel('$x_2$','interpreter','latex','fontsize',fs)
    title('$ \| \textbf{v}_{(0)} \|_2 $','interpreter','latex','fontsize',fs)
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize=0.6*fs;
    
    subplot(312)
    contourf(x,y,sqrt(P1.o1.up(:,:,j)'.^2+P1.o1.vp(:,:,j)'.^2),levels,'edgecolor','none'); hold on
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
    xlabel('$x_1$','interpreter','latex','fontsize',fs)
    ylabel('$x_2$','interpreter','latex','fontsize',fs)
    title('$ \| \textbf{v}_{(1)} \|_2 $','interpreter','latex','fontsize',fs)
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize=0.6*fs;
    
    
    subplot(313)
    contourf(x,y,sqrt(P1.o2.up(:,:,j)'.^2+P1.o2.vp(:,:,j)'.^2),levels,'edgecolor','none'); hold on
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
    xlabel('$x_1$','interpreter','latex','fontsize',fs)
    ylabel('$x_2$','interpreter','latex','fontsize',fs)
    title('$ \| \textbf{v}_{(2)} \|_2 $','interpreter','latex','fontsize',fs)
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize=0.6*fs;
    
    
    drawnow
    
    j
end


%% PREF CONCENTRATION =====================================================

close all
figure('units','inches','outerposition',[1 1 8+1 8+1],'color','w')

for j=1:1:nt(1)
    
    subplot(321)
    xx = P1.o0.xp(:,:,j);
    yy = P1.o0.yp(:,:,j);
    scatter(xx(:),yy(:),0.01,'k')
    axis equal
    hold off
%     xlim([0 2])
%     ylim([0 1])
    ax=gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize=0.6*fs;
%     ax.XTick = [];
%     ax.YTick = [];
    % title('$St = 0.5, \ \sigma_{f}=0$','interpreter','latex','fontsize',fs)
    
    subplot(323)
    xx = P1.o1.xp(:,:,j);
    yy = P1.o1.yp(:,:,j);
    scatter(xx(:),yy(:),0.01,'k')
    axis equal
    hold off
%     xlim([0 2])
%     ylim([0 1])
    ax=gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize=0.6*fs;
%     ax.XTick = [];
%     ax.YTick = [];
    % title('$St = 0.5, \ \sigma_{f}=0$','interpreter','latex','fontsize',fs)
    
    subplot(325)
    xx = P1.o2.xp(:,:,j);
    yy = P1.o2.yp(:,:,j);
    scatter(xx(:),yy(:),0.01,'k')
    axis equal
    hold off
%     xlim([0 2])
%     ylim([0 1])
    ax=gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize=0.6*fs;
%     ax.XTick = [];
%     ax.YTick = [];
    % title('$St = 0.5, \ \sigma_{f}=0$','interpreter','latex','fontsize',fs)

    subplot(322)
    xx = P1.o0.up(:,:,j);
    yy = P1.o0.vp(:,:,j);
    scatter(xx(:),yy(:),0.01,'k')
    axis equal
    hold off
%     xlim([0 2])
%     ylim([0 1])
    ax=gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize=0.6*fs;
%     ax.XTick = [];
%     ax.YTick = [];
    % title('$St = 0.5, \ \sigma_{f}=0$','interpreter','latex','fontsize',fs)
    
    subplot(324)
    xx = P1.o1.up(:,:,j);
    yy = P1.o1.vp(:,:,j);
    scatter(xx(:),yy(:),0.01,'k')
    axis equal
    hold off
%     xlim([0 2])
%     ylim([0 1])
    ax=gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize=0.6*fs;
%     ax.XTick = [];
%     ax.YTick = [];
    % title('$St = 0.5, \ \sigma_{f}=0$','interpreter','latex','fontsize',fs)
    
    subplot(326)
    xx = P1.o2.up(:,:,j);
    yy = P1.o2.vp(:,:,j);
    scatter(xx(:),yy(:),0.01,'k')
    axis equal
    hold off
%     xlim([0 2])
%     ylim([0 1])
    ax=gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize=0.6*fs;
%     ax.XTick = [];
%     ax.YTick = [];
    % title('$St = 0.5, \ \sigma_{f}=0$','interpreter','latex','fontsize',fs)

    
    
    drawnow
end

% if saveFigs == 1
%     cd([figuresDir '_' label]);
%     thisName = 'dG_intro_part_St1_a0_T12w6';
%     %     savefig(iFig,thisName)
%     %     print(iFig,thisName,'-depsc','-painters')
%     print(iFig,thisName,'-dpng','-r600','-painters')
%     cd ..;
% end
% iFig = iFig+1; clear p1 leg1 leg







%%
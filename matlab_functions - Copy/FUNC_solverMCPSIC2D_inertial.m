%% ========================================================================

function P1 = FUNC_solverMCPSIC2D_inertial(xp0,yp0,up0,vp0,tLim,nt, ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup)

disp('FUNC_solverMCPSIC2D_inertial')

% Small calculations ------------------------------------------------------
nt_2save = nt(1);
nt_2iter = nt(2);
Nt       = (nt(1)-1)*(nt(2)-1)+1;
dt       = (tLim(2)-tLim(1))/(Nt-1);
t        = linspace(tLim(1),tLim(2),nt_2save)';

% Initializating ----------------------------------------------------------
npx  = length(xp0(:,1,1)); % clouds along x
npy  = length(xp0(1,:,1)); % clouds along y
ns   = length(xp0(1,1,:)); % samples per cloud
np   = npx*npy*ns; % total clouds
xp0c = xp0(:); % grid of initial particles put into columns
yp0c = yp0(:);
up0c = up0(:);
vp0c = vp0(:);
time = tLim(1);
yy   = [ xp0c ; yp0c ; up0c ; vp0c ];
y    = yy(:,1);
% RHS(:,1) = FUNC_rhsMCPSIC2D_inertial(time,y,flowType,dataFlow,f1Type,dataf1,npx,npy,ns,taup);

% Loop in time ------------------------------------------------------------
for j=1:nt_2save-1            % loop used to save in time
    for dummyVar=1:nt_2iter-1 % loop that iterates in every time step to save
        switch timeMethod
            case 1 % EE
                k1 = FUNC_rhsMCPSIC2D_inertial(time       ,                        y,flowType,dataFlow,f1Type,dataf1,npx,npy,ns,taup);
                y  = y+dt*k1;
            case 2 % RK2
                disp('not coded')
            case 3 % RK3 TVD
                k1 = FUNC_rhsMCPSIC2D_inertial(time       ,                        y,flowType,dataFlow,f1Type,dataf1,npx,npy,ns,taup);
                k2 = FUNC_rhsMCPSIC2D_inertial(time+dt    ,                  y+dt*k1,flowType,dataFlow,f1Type,dataf1,npx,npy,ns,taup);
                k3 = FUNC_rhsMCPSIC2D_inertial(time+3/4*dt,y+(1/4)*dt*k1+(1/4)*dt*k2,flowType,dataFlow,f1Type,dataf1,npx,npy,ns,taup);
                y  = y+dt*((1/6)*k1+(1/6)*k2+(2/3)*k3);
            case 4 % RK4
                k1 = FUNC_rhsMCPSIC2D_inertial(time       ,y           ,flowType,dataFlow,f1Type,dataf1,npx,npy,ns,taup);
                k2 = FUNC_rhsMCPSIC2D_inertial(time+0.5*dt,y+0.5*dt.*k1,flowType,dataFlow,f1Type,dataf1,npx,npy,ns,taup);
                k3 = FUNC_rhsMCPSIC2D_inertial(time+0.5*dt,y+0.5*dt.*k2,flowType,dataFlow,f1Type,dataf1,npx,npy,ns,taup);
                k4 = FUNC_rhsMCPSIC2D_inertial(time+dt    ,y+    dt.*k3,flowType,dataFlow,f1Type,dataf1,npx,npy,ns,taup);
                y  = y+dt*(k1+2*k2+2*k3+k4)/6;
            otherwise
                disp('Error in timeMethod')
        end
        time = time + dt;
    end
    yy(:,j+1) = y; % Saving after iterating without saving
%     RHS(:,j+1) = k1;
    if j==round(0.1*nt(1))
        fprintf('10%%..  t=%1.4f \n',time);
    elseif j==round(0.25*nt(1))
        fprintf('25%%..  t=%1.4f \n',time);
    elseif j==round(0.5*nt(1))
        fprintf('50%%..  t=%1.4f \n',time);
    elseif j==round(0.75*nt(1))
        fprintf('75%%..  t=%1.4f \n',time);
    elseif j==round(nt(1)-1)
        fprintf('100%%.  t=%1.4f \n',time);
    end
end

% Reordering and extracting -----------------------------------------------
P1.xp = reshape(yy(0*np+1:1*np,:),npx,npy,ns,nt(1));
P1.yp = reshape(yy(1*np+1:2*np,:),npx,npy,ns,nt(1));
P1.up = reshape(yy(2*np+1:3*np,:),npx,npy,ns,nt(1));
P1.vp = reshape(yy(3*np+1:4*np,:),npx,npy,ns,nt(1));
P1.t  = t;

% P1.RHS.xp = reshape(RHS(0*np+1:1*np,:),npx,npy,ns,nt(1));
% P1.RHS.yp = reshape(RHS(1*np+1:2*np,:),npx,npy,ns,nt(1));

% Post --------------------------------------------------------------------
for j=1:nt(1)
    for ix=1:npx
        for iy=1:npy
            P1.mean_xp(ix,iy,j)   = mean(squeeze(P1.xp(ix,iy,:,j)));
            P1.mean_yp(ix,iy,j)   = mean(squeeze(P1.yp(ix,iy,:,j)));
            P1.mean_up(ix,iy,j)   = mean(squeeze(P1.up(ix,iy,:,j)));
            P1.mean_vp(ix,iy,j)   = mean(squeeze(P1.vp(ix,iy,:,j)));
            P1.xpxp(ix,iy,j)      = mean((squeeze(P1.xp(ix,iy,:,j))-P1.mean_xp(ix,iy,j)).^2);
            P1.xpyp(ix,iy,j)      = mean((squeeze(P1.xp(ix,iy,:,j))-P1.mean_xp(ix,iy,j)).*(squeeze(P1.yp(ix,iy,:,j))-P1.mean_yp(ix,iy,j)));
            P1.xpup(ix,iy,j)      = mean((squeeze(P1.xp(ix,iy,:,j))-P1.mean_xp(ix,iy,j)).*(squeeze(P1.up(ix,iy,:,j))-P1.mean_up(ix,iy,j)));
            P1.xpvp(ix,iy,j)      = mean((squeeze(P1.xp(ix,iy,:,j))-P1.mean_xp(ix,iy,j)).*(squeeze(P1.vp(ix,iy,:,j))-P1.mean_vp(ix,iy,j)));
            P1.ypyp(ix,iy,j)      = mean((squeeze(P1.yp(ix,iy,:,j))-P1.mean_yp(ix,iy,j)).^2);
            P1.ypup(ix,iy,j)      = mean((squeeze(P1.yp(ix,iy,:,j))-P1.mean_yp(ix,iy,j)).*(squeeze(P1.up(ix,iy,:,j))-P1.mean_up(ix,iy,j)));
            P1.ypvp(ix,iy,j)      = mean((squeeze(P1.yp(ix,iy,:,j))-P1.mean_yp(ix,iy,j)).*(squeeze(P1.vp(ix,iy,:,j))-P1.mean_vp(ix,iy,j)));
            P1.upup(ix,iy,j)      = mean((squeeze(P1.up(ix,iy,:,j))-P1.mean_up(ix,iy,j)).^2);
            P1.upvp(ix,iy,j)      = mean((squeeze(P1.up(ix,iy,:,j))-P1.mean_up(ix,iy,j)).*(squeeze(P1.vp(ix,iy,:,j))-P1.mean_vp(ix,iy,j)));
            P1.vpvp(ix,iy,j)      = mean((squeeze(P1.vp(ix,iy,:,j))-P1.mean_vp(ix,iy,j)).^2);
            
            % Third
            P1.xpxpxp(ix,iy,j)    = mean((squeeze(P1.xp(ix,iy,:,j))-P1.mean_xp(ix,iy,j)).^3);
            P1.upupup(ix,iy,j)    = mean((squeeze(P1.up(ix,iy,:,j))-P1.mean_up(ix,iy,j)).^3);
                        
            % Positions
            Kxp                   = [P1.xpxp(ix,iy,j) P1.xpyp(ix,iy,j) ; P1.xpyp(ix,iy,j) P1.ypyp(ix,iy,j)];
            if sum(isnan(Kxp(:))) > 0 || sum(isinf(Kxp(:))) > 0
                P1.maxEigKxp(ix,iy,j) = 0;
            else
                P1.maxEigKxp(ix,iy,j) = max(eig(Kxp));
            end
            T                     = abs(P1.t(j)-tLim(1));
            P1.FTLE(ix,iy,j)      = 1/T*log(sqrt(P1.maxEigKxp(ix,iy,j)/P1.maxEigKxp(ix,iy,1)));
            % Velocities
            Kup                   = [P1.upup(ix,iy,j) P1.upvp(ix,iy,j) ; P1.upvp(ix,iy,j) P1.vpvp(ix,iy,j)];
            if sum(isnan(Kup(:))) > 0 || sum(isinf(Kup(:))) > 0
                P1.maxEigKup(ix,iy,j) = 0;
            else
                P1.maxEigKup(ix,iy,j) = max(eig(Kup));
            end
%             P1.FTLEv(ix,iy,j)     = (1/(P1.t(j)-tLim(1)))*log(sqrt(P1.maxEigKup(ix,iy,j)/P1.maxEigKup(ix,iy,1)));
            P1.FTLEv(ix,iy,j)     = 1/T*log(sqrt(P1.maxEigKup(ix,iy,j)));
            % Erros
            P1.mu1(ix,iy,j)       = sqrt(P1.mean_xp(ix,iy,j).^2+P1.mean_yp(ix,iy,j).^2+P1.mean_up(ix,iy,j).^2+P1.mean_vp(ix,iy,j).^2);
            K                     = [ ...
                P1.xpxp(ix,iy,j) P1.xpyp(ix,iy,j) P1.xpup(ix,iy,j) P1.xpvp(ix,iy,j) ;
                P1.xpyp(ix,iy,j) P1.ypyp(ix,iy,j) P1.ypup(ix,iy,j) P1.ypvp(ix,iy,j) ; 
                P1.xpup(ix,iy,j) P1.ypup(ix,iy,j) P1.upup(ix,iy,j) P1.upvp(ix,iy,j) ; 
                P1.xpvp(ix,iy,j) P1.ypvp(ix,iy,j) P1.upvp(ix,iy,j) P1.vpvp(ix,iy,j)];
            P1.mu2(ix,iy,j)       = det(K);
%             % Relative velocity magnitudes
%             out1 = FUNC_flow(time,squeeze(P1.xp(ix,iy,:,j)),squeeze(P1.yp(ix,iy,:,j)),squeeze(P1.xp(ix,iy,:,j))*0,flowType,dataFlow);
%             u    = out1.u;
%             v    = out1.v;
%             ax   = u - squeeze(P1.up(ix,iy,:,j));
%             ay   = v - squeeze(P1.vp(ix,iy,:,j));
%             axax = mean((ax-mean(ax)).^2);
%             axay = mean((ax-mean(ax)).*(ay-mean(ay)));
%             ayay = mean((ay-mean(ay)).^2);
%             Ka   = [axax axay ; axay ayay];
%             if sum(isnan(Ka(:))) > 0 || sum(isinf(Ka(:))) > 0
%                 P1.maxEigKa(ix,iy,j) = 0;
%             else
%                 P1.maxEigKa(ix,iy,j) = max(eig(Ka));
%             end
        end
    end
%     j
end
P1.FTLE(:,:,1)  = 0;
P1.FTLEv(:,:,1) = 0;

% for j=1:nt(1)
%     out.cx(:,j)  = [ out.xp(:,j) ; out.xp(1,j) ];
%     out.cy(:,j)  = [ out.yp(:,j) ; out.yp(1,j) ];
%     [out.curv_L(:,j),out.curv_R(:,j),out.curv_K(:,:,j)] = curvature([out.cx(:,j) out.cy(:,j)]);
%     [out.curv_maxR(j),out.curv_maxR_ind(j)]             = max(1./out.curv_R(:,j));
% end




end









%% ========================================================================

function S1 = FUNC_solverSPARSER2D_inertial(mean_xp0,mean_yp0,mean_up0,mean_vp0, ...
    xpxp0,xpyp0,xpup0,xpvp0,ypyp0,ypup0,ypvp0,upup0,upvp0,vpvp0, ...
    tLim,nt,flowType,dataFlow,f1Type,dataf1,timeMethod,nmps,taup)

disp('FUNC_solverSPARSE2D_inertial')

% Small calculations ------------------------------------------------------
nt_2save = nt(1);
nt_2iter = nt(2);
Nt       = (nt(1)-1)*(nt(2)-1)+1;
dt       = (tLim(2)-tLim(1))/(Nt-1);
t        = linspace(tLim(1),tLim(2),nt_2save)';
npx      = length(mean_xp0(:,1)); % clouds along x
npy      = length(mean_xp0(1,:)); % clouds along y
np       = npx*npy; % total clouds

% Splitting ---------------------------------------------------------------
nmp = nmps(1)*nmps(2);
M   = max([nmps(1) nmps(2)]); % The splitting level is constant
NP  = np*nmp;
fprintf('Mp per cloud: %1.0f, total split solution: %1.0f \n',nmp,NP)

% Prepare the parameter matrices for each macro-particle to the
% SPARSE simulation
for jx=1:npx
    for jy=1:npy
        k = 1;
        for ix=1:nmps(1)
            for iy=1:nmps(2)
                if M == 1
                    Lxp = 0; Dxp = 0;
                    Lyp = 0; Dyp = 0;
                elseif M == 2
                    Lxp = 2*sqrt( 0.5*xpxp0(jx,jy) ); Dxp = 0.5*Lxp;
                    Lyp = 2*sqrt( 0.5*ypyp0(jx,jy) ); Dyp = 0.5*Lyp;
                elseif M == 3
                    Lxp = sqrt( xpxp0(jx,jy) ); Dxp = Lxp;
                    Lyp = sqrt( ypyp0(jx,jy) ); Dyp = Lyp;
                elseif M == 4
                    Lxp = sqrt( 3*xpxp0(jx,jy)/(2*1.5^2+2*0.5^2) ); Dxp = 1.5*Lxp;
                    Lyp = sqrt( 3*ypyp0(jx,jy)/(2*1.5^2+2*0.5^2) ); Dyp = 1.5*Lyp;
                elseif M == 5
                    Lxp = sqrt( 4/10 *xpxp0(jx,jy) ); Dxp = 2*Lxp;
                    Lyp = sqrt( 4/10 *ypyp0(jx,jy) ); Dyp = 2*Lyp;
                elseif M == 6
                    Lxp = sqrt( 5/2 *xpxp0(jx,jy)/(2.5^2+1.5^2+0.5^2) ); Dxp = 2.5*Lxp;
                    Lyp = sqrt( 5/2 *ypyp0(jx,jy)/(2.5^2+1.5^2+0.5^2) ); Dyp = 2.5*Lyp;
                elseif M == 7
                    Lxp = sqrt( 6/2 *xpxp0(jx,jy)/(3^2+2^2+1^2) ); Dxp = 3*Lxp;
                    Lyp = sqrt( 6/2 *ypyp0(jx,jy)/(3^2+2^2+1^2) ); Dyp = 3*Lyp;
                elseif M == 8
                    Lxp = sqrt( 7/2 *xpxp0(jx,jy)/(3.5^2+2.5^2+1.5^2+0.5^2) ); Dxp = 3.5*Lxp;
                    Lyp = sqrt( 7/2 *ypyp0(jx,jy)/(3.5^2+2.5^2+1.5^2+0.5^2) ); Dyp = 3.5*Lyp;
                end
                mean_xp0s(jx,jy,k) = mean_xp0(jx,jy) - Dxp + (ix-1)*Lxp;
                mean_yp0s(jx,jy,k) = mean_yp0(jx,jy) - Dyp + (iy-1)*Lyp;
                xpxp0s(jx,jy,k)    = xpxp0(jx,jy)/M;
                ypyp0s(jx,jy,k)    = ypyp0(jx,jy)/M;
                xpyp0s(jx,jy,k)    = xpyp0(jx,jy);
                % Along u and v ther eis no splitting, assumed zero
                mean_up0s(jx,jy,k) = 0;
                mean_vp0s(jx,jy,k) = 0;
                xpup0s(jx,jy,k)    = 0;
                xpvp0s(jx,jy,k)    = 0;
                ypup0s(jx,jy,k)    = 0;
                ypvp0s(jx,jy,k)    = 0;
                upup0s(jx,jy,k)    = 0;
                upvp0s(jx,jy,k)    = 0;
                vpvp0s(jx,jy,k)    = 0;
                
                k = k + 1;
            end
        end
    end
end

% Initializating ----------------------------------------------------------
mean_xp0c = mean_xp0s(:); % grid of initial clouds put into columns
mean_yp0c = mean_yp0s(:);
mean_up0c = mean_up0s(:);
mean_vp0c = mean_vp0s(:);
xpxp0c    = xpxp0s(:);
xpyp0c    = xpyp0s(:);
xpup0c    = xpup0s(:);
xpvp0c    = xpvp0s(:);
ypyp0c    = ypyp0s(:);
ypup0c    = ypup0s(:);
ypvp0c    = ypvp0s(:);
upup0c    = upup0s(:);
upvp0c    = upvp0s(:);
vpvp0c    = vpvp0s(:);
time      = tLim(1);
yy        = [ mean_xp0c ; mean_yp0c ; mean_up0c ; mean_vp0c ; ...
    xpxp0c ; xpyp0c ; xpup0c ; xpvp0c ; ...
    ypyp0c ; ypup0c ; ypvp0c ; ...
    upup0c ; upvp0c ; ...
    vpvp0c ];
y         = yy(:,1);

% Loop in time ------------------------------------------------------------
for j=1:nt_2save-1            % loop used to save in time
    for dummyVar=1:nt_2iter-1 % loop that iterates in every time step to save
        switch timeMethod
            case 1 % EE
                k1 = FUNC_rhsSPARSE2D_inertial(time,y,flowType,dataFlow,f1Type,dataf1,npx,npy,nmps,taup);
                y  = y+dt*k1;
            case 2 % RK2
                disp('RK2 not coded')
            case 3 % RK3 TVD
                k1 = FUNC_rhsSPARSE2D_inertial(time       ,                        y,flowType,dataFlow,f1Type,dataf1,npx,npy,nmps,taup);
                k2 = FUNC_rhsSPARSE2D_inertial(time+dt    ,                  y+dt*k1,flowType,dataFlow,f1Type,dataf1,npx,npy,nmps,taup);
                k3 = FUNC_rhsSPARSE2D_inertial(time+3/4*dt,y+(1/4)*dt*k1+(1/4)*dt*k2,flowType,dataFlow,f1Type,dataf1,npx,npy,nmps,taup);
                y  = y+dt*((1/6)*k1+(1/6)*k2+(2/3)*k3);
            case 4 % RK4
                k1 = FUNC_rhsSPARSE2D_inertial(time       ,y           ,flowType,dataFlow,f1Type,dataf1,npx,npy,nmps,taup);
                k2 = FUNC_rhsSPARSE2D_inertial(time+0.5*dt,y+0.5*dt.*k1,flowType,dataFlow,f1Type,dataf1,npx,npy,nmps,taup);
                k3 = FUNC_rhsSPARSE2D_inertial(time+0.5*dt,y+0.5*dt.*k2,flowType,dataFlow,f1Type,dataf1,npx,npy,nmps,taup);
                k4 = FUNC_rhsSPARSE2D_inertial(time+dt    ,y+    dt.*k3,flowType,dataFlow,f1Type,dataf1,npx,npy,nmps,taup);
                y  = y+dt*(k1+2*k2+2*k3+k4)/6;
            otherwise
                disp('Error in timeMethod')
        end
        time = time + dt;
    end
    yy(:,j+1)  = y; % Saving after iterating without saving
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
mean_xp = reshape(yy( 0*NP+1: 1*NP,:),npx,npy,nmp,nt(1));
mean_yp = reshape(yy( 1*NP+1: 2*NP,:),npx,npy,nmp,nt(1));
mean_up = reshape(yy( 2*NP+1: 3*NP,:),npx,npy,nmp,nt(1));
mean_vp = reshape(yy( 3*NP+1: 4*NP,:),npx,npy,nmp,nt(1));
xpxp    = reshape(yy( 4*NP+1: 5*NP,:),npx,npy,nmp,nt(1));
xpyp    = reshape(yy( 5*NP+1: 6*NP,:),npx,npy,nmp,nt(1));
xpup    = reshape(yy( 6*NP+1: 7*NP,:),npx,npy,nmp,nt(1));
xpvp    = reshape(yy( 7*NP+1: 8*NP,:),npx,npy,nmp,nt(1));
ypyp    = reshape(yy( 8*NP+1: 9*NP,:),npx,npy,nmp,nt(1));
ypup    = reshape(yy( 9*NP+1:10*NP,:),npx,npy,nmp,nt(1));
ypvp    = reshape(yy(10*NP+1:11*NP,:),npx,npy,nmp,nt(1));
upup    = reshape(yy(11*NP+1:12*NP,:),npx,npy,nmp,nt(1));
upvp    = reshape(yy(12*NP+1:13*NP,:),npx,npy,nmp,nt(1));
vpvp    = reshape(yy(13*NP+1:14*NP,:),npx,npy,nmp,nt(1));
S1.t       = t;

% Joining macro-particles -------------------------------------------------
w = 1/nmp; % the weights would be all equal
for k=1:nt(1)
    for i=1:npx
        for j=1:npy
            smxp  = squeeze(mean_xp(i,j,:,k));
            smyp  = squeeze(mean_yp(i,j,:,k));
            smup  = squeeze(mean_up(i,j,:,k));
            smvp  = squeeze(mean_vp(i,j,:,k));
            sxpxp = squeeze(xpxp(i,j,:,k));
            sxpyp = squeeze(xpyp(i,j,:,k));
            sxpup = squeeze(xpup(i,j,:,k));
            sxpvp = squeeze(xpvp(i,j,:,k));
            sypyp = squeeze(ypyp(i,j,:,k));
            sypup = squeeze(ypup(i,j,:,k));
            sypvp = squeeze(ypvp(i,j,:,k));
            supup = squeeze(upup(i,j,:,k));
            supvp = squeeze(upvp(i,j,:,k));
            svpvp = squeeze(vpvp(i,j,:,k));
            S1.mean_xp(i,j,k) = sum(w*smxp);
            S1.mean_yp(i,j,k) = sum(w*smyp);
            S1.mean_up(i,j,k) = sum(w*smup);
            S1.mean_vp(i,j,k) = sum(w*smvp);
            S1.xpxp(i,j,k)    = sum(w*sxpxp) + sum(w*(smxp - S1.mean_xp(i,j,k)).^2);
            S1.xpyp(i,j,k)    = sum(w*sxpyp) + sum(w*(smxp - S1.mean_xp(i,j,k)).*(smyp - S1.mean_yp(i,j,k)));
            S1.xpup(i,j,k)    = sum(w*sxpup) + sum(w*(smxp - S1.mean_xp(i,j,k)).*(smup - S1.mean_up(i,j,k)));
            S1.xpvp(i,j,k)    = sum(w*sxpvp) + sum(w*(smxp - S1.mean_xp(i,j,k)).*(smvp - S1.mean_vp(i,j,k)));
            S1.ypyp(i,j,k)    = sum(w*sypyp) + sum(w*(smyp - S1.mean_yp(i,j,k)).^2);
            S1.ypup(i,j,k)    = sum(w*sypup) + sum(w*(smyp - S1.mean_yp(i,j,k)).*(smup - S1.mean_up(i,j,k)));
            S1.ypvp(i,j,k)    = sum(w*sypvp) + sum(w*(smyp - S1.mean_yp(i,j,k)).*(smvp - S1.mean_vp(i,j,k)));
            S1.upup(i,j,k)    = sum(w*supup) + sum(w*(smup - S1.mean_up(i,j,k)).^2);
            S1.upvp(i,j,k)    = sum(w*supvp) + sum(w*(smup - S1.mean_up(i,j,k)).*(smvp - S1.mean_vp(i,j,k)));
            S1.vpvp(i,j,k)    = sum(w*svpvp) + sum(w*(smvp - S1.mean_vp(i,j,k)).^2);
        end
    end
end

% Post --------------------------------------------------------------------
for j=1:nt(1)
    for ix=1:npx
        for iy=1:npy
            Kxp                   = [S1.xpxp(ix,iy,j) S1.xpyp(ix,iy,j) ; S1.xpyp(ix,iy,j) S1.ypyp(ix,iy,j)];
            S1.maxEigKxp(ix,iy,j) = max(eig(Kxp));
            S1.FTLE(ix,iy,j)      = (1/(S1.t(j)-tLim(1)))*log(sqrt(S1.maxEigKxp(ix,iy,j)/S1.maxEigKxp(ix,iy,1)));
            Kup                   = [S1.upup(ix,iy,j) S1.upvp(ix,iy,j) ; S1.upvp(ix,iy,j) S1.vpvp(ix,iy,j)];
            S1.maxEigKup(ix,iy,j) = max(eig(Kup));
%             S1.FTLEv(ix,iy,j)     = (1/(S1.t(j)-tLim(1)))*log(sqrt(S1.maxEigKup(ix,iy,j)/S1.maxEigKup(ix,iy,1)));
            S1.FTLEv(ix,iy,j)     = (1/(S1.t(j)-tLim(1)))*log(sqrt(S1.maxEigKup(ix,iy,j)));
            S1.mu1(ix,iy,j)       = sqrt(S1.mean_xp(ix,iy,j).^2+S1.mean_yp(ix,iy,j).^2+S1.mean_up(ix,iy,j).^2+S1.mean_vp(ix,iy,j).^2);
            K                     = [ ...
                S1.xpxp(ix,iy,j) S1.xpyp(ix,iy,j) S1.xpup(ix,iy,j) S1.xpvp(ix,iy,j) ;
                S1.xpyp(ix,iy,j) S1.ypyp(ix,iy,j) S1.ypup(ix,iy,j) S1.ypvp(ix,iy,j) ;
                S1.xpup(ix,iy,j) S1.ypup(ix,iy,j) S1.upup(ix,iy,j) S1.upvp(ix,iy,j) ;
                S1.xpvp(ix,iy,j) S1.ypvp(ix,iy,j) S1.upvp(ix,iy,j) S1.vpvp(ix,iy,j)];
            S1.mu2(ix,iy,j)       = det(K);
        end
    end
end
S1.FTLE(:,:,1)  = 0;
S1.FTLEv(:,:,1) = 0;
S1.FTLE         = real(S1.FTLE);
S1.FTLEv        = real(S1.FTLEv);

end
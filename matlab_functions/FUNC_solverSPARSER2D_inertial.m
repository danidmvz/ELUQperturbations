%% ========================================================================

function S1 = FUNC_solverSPARSER2D_inertial(mean_xp0,mean_yp0,mean_up0,mean_vp0, ...
    mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,M,taup)

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
NP  = np*M;
nmp = M;
fprintf('nmp per cloud: %1.0f, total split solution: %1.0f \n',M,NP)

% Prepare the parameter matrices for each macro-particle to the
% SPARSE simulation
for ix=1:npx
    for iy=1:npy
        k = 1;
        for ia=1:M
            if M == 1
                La = 0;
                Da = 0;
            elseif M == 2
                La = 2*sqrt( 0.5*aa );
                Da = 0.5*La;
            elseif M == 3
                La = sqrt( aa );
                Da = La;
            elseif M == 4
                La = sqrt( 3*aa/(2*1.5^2+2*0.5^2) );
                Da = 1.5*La;
            elseif M == 5
                La = sqrt( 4/10 *aa );
                Da = 2*La;
            elseif M == 6
                La = sqrt( 5/2 *aa/(2.5^2+1.5^2+0.5^2) );
                Da = 2.5*La;
            elseif M == 7
                La = sqrt( 6/2 *aa/(3^2+2^2+1^2) );
                Da = 3*La;
            elseif M == 8
                La = sqrt( 7/2 *aa/(3.5^2+2.5^2+1.5^2+0.5^2) );
                Da = 3.5*La;
            end
            mean_a_s(ix,iy,k)  = mean_a - Da + (ia-1)*La;
            aa_s(ix,iy,k)      = aa/M;
            mean_xp0s(ix,iy,k) = mean_xp0(ix,iy);
            mean_yp0s(ix,iy,k) = mean_yp0(ix,iy);
            mean_up0s(ix,iy,k) = mean_up0(ix,iy);
            mean_vp0s(ix,iy,k) = mean_vp0(ix,iy);
            xpxp0s(ix,iy,k)    = 0;
            xpyp0s(ix,iy,k)    = 0;
            xpup0s(ix,iy,k)    = 0;
            xpvp0s(ix,iy,k)    = 0;
            ypyp0s(ix,iy,k)    = 0;
            ypup0s(ix,iy,k)    = 0;
            ypvp0s(ix,iy,k)    = 0;
            upup0s(ix,iy,k)    = 0;
            upvp0s(ix,iy,k)    = 0;
            vpvp0s(ix,iy,k)    = 0;
            axp_s(ix,iy,k)     = 0;
            axp_s(ix,iy,k)     = 0;
            ayp_s(ix,iy,k)     = 0;
            aup_s(ix,iy,k)     = 0;
            avp_s(ix,iy,k)     = 0;
            k = k + 1;
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
axp0c     = axp_s(:);
ayp0c     = ayp_s(:);
aup0c     = aup_s(:);
avp0c     = avp_s(:);
time      = tLim(1);
yy        = [ mean_xp0c ; mean_yp0c ; mean_up0c ; mean_vp0c ; ...
    xpxp0c ; xpyp0c ; xpup0c ; xpvp0c ; ...
    ypyp0c ; ypup0c ; ypvp0c ; ...
    upup0c ; upvp0c ; ...
    vpvp0c ; 
    axp0c  ; ayp0c  ; aup0c  ; avp0c];
y         = yy(:,1);

% Loop in time ------------------------------------------------------------
for j=1:nt_2save-1            % loop used to save in time
    for dummyVar=1:nt_2iter-1 % loop that iterates in every time step to save
        switch timeMethod
            case 1 % EE
                k1 = FUNC_rhsSPARSER2D_inertial(time,y,flowType,dataFlow,taup,mean_a_s,aa_s,npx,npy,M);
                y  = y+dt*k1;     

            case 2 % RK2
                disp('RK2 not coded')
            case 3 % RK3 TVD
                k1 = FUNC_rhsSPARSER2D_inertial(time       ,                        y,flowType,dataFlow,taup,mean_a_s,aa_s,npx,npy,M);
                k2 = FUNC_rhsSPARSER2D_inertial(time+dt    ,                  y+dt*k1,flowType,dataFlow,taup,mean_a_s,aa_s,npx,npy,M);
                k3 = FUNC_rhsSPARSER2D_inertial(time+3/4*dt,y+(1/4)*dt*k1+(1/4)*dt*k2,flowType,dataFlow,taup,mean_a_s,aa_s,npx,npy,M);
                y  = y+dt*((1/6)*k1+(1/6)*k2+(2/3)*k3);
            case 4 % RK4
                k1 = FUNC_rhsSPARSER2D_inertial(time       ,y           ,flowType,dataFlow,taup,mean_a_s,aa_s,npx,npy,M);
                k2 = FUNC_rhsSPARSER2D_inertial(time+0.5*dt,y+0.5*dt.*k1,flowType,dataFlow,taup,mean_a_s,aa_s,npx,npy,M);
                k3 = FUNC_rhsSPARSER2D_inertial(time+0.5*dt,y+0.5*dt.*k2,flowType,dataFlow,taup,mean_a_s,aa_s,npx,npy,M);
                k4 = FUNC_rhsSPARSER2D_inertial(time+dt    ,y+    dt.*k3,flowType,dataFlow,taup,mean_a_s,aa_s,npx,npy,M);
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
axp     = reshape(yy(14*NP+1:15*NP,:),npx,npy,nmp,nt(1));
ayp     = reshape(yy(15*NP+1:16*NP,:),npx,npy,nmp,nt(1));
aup     = reshape(yy(16*NP+1:17*NP,:),npx,npy,nmp,nt(1));
avp     = reshape(yy(17*NP+1:18*NP,:),npx,npy,nmp,nt(1));
S1.t    = t;

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
            saxp  = squeeze(axp(i,j,:,k));
            sayp  = squeeze(ayp(i,j,:,k));
            saup  = squeeze(aup(i,j,:,k));
            savp  = squeeze(avp(i,j,:,k));
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
            S1.axp(i,j,k)     = sum(w*saxp) + sum(w*(squeeze(mean_a_s(i,j,:)) - mean_a).*(smxp - S1.mean_xp(i,j,k)));
            S1.ayp(i,j,k)     = sum(w*sayp) + sum(w*(squeeze(mean_a_s(i,j,:)) - mean_a).*(smyp - S1.mean_yp(i,j,k)));
            S1.aup(i,j,k)     = sum(w*saup) + sum(w*(squeeze(mean_a_s(i,j,:)) - mean_a).*(smup - S1.mean_up(i,j,k)));
            S1.avp(i,j,k)     = sum(w*savp) + sum(w*(squeeze(mean_a_s(i,j,:)) - mean_a).*(smvp - S1.mean_vp(i,j,k)));
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
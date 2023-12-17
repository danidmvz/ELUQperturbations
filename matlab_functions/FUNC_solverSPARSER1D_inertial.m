%% ========================================================================

function S1 = FUNC_solverSPARSER1D_inertial(mean_xp0,mean_up0, ...
    mean_a,aa,tLim,nt,flowType,dataFlow,timeMethod,M,taup)

disp('FUNC_solverSPARSER1D_inertial')

% Small calculations ------------------------------------------------------
nt_2save = nt(1);
nt_2iter = nt(2);
Nt       = (nt(1)-1)*(nt(2)-1)+1;
dt       = (tLim(2)-tLim(1))/(Nt-1);
t        = linspace(tLim(1),tLim(2),nt_2save)';

% Splitting ---------------------------------------------------------------
nmp = M; % There are 7 moments per cloud btu I only split along \alpha
fprintf('nmp per cloud: %1.0f \n',nmp)

% Prepare the parameter matrices for each macro-particle to the
% SPARSE simulation
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
    mean_a_s(k,1)  = mean_a - Da + (ia-1)*La;
    aa_s(k,1)      = aa/M;
    mean_xp0s(k,1) = mean_xp0;
    mean_up0s(k,1) = mean_up0;
    xpxp0s(k,1)    = 0;
    xpup0s(k,1)    = 0;
    upup0s(k,1)    = 0;
    axp_s(k,1)     = 0;
    aup_s(k,1)     = 0;
    k = k + 1;
end

% Initializating ----------------------------------------------------------
time = tLim(1);
yy   = [ mean_xp0s(:) ; mean_up0s(:) ; xpxp0s(:) ; xpup0s(:) ; upup0s(:) ; ...
    axp_s(:); aup_s(:)];
y    = yy(:,1);

% Loop in time ------------------------------------------------------------
for j=1:nt_2save-1            % loop used to save in time
    for dummyVar=1:nt_2iter-1 % loop that iterates in every time step to save
        switch timeMethod
            case 1 % EE
                k1 = FUNC_rhsSPARSER1D_inertial(time,y,flowType,dataFlow,taup,mean_a_s,aa_s);
                y  = y+dt*k1;
            case 2 % RK2
                disp('RK2 not coded')
            case 3 % RK3 TVD
                k1 = FUNC_rhsSPARSER1D_inertial(time       ,                        y,flowType,dataFlow,taup,mean_a_s,aa_s);
                k2 = FUNC_rhsSPARSER1D_inertial(time+dt    ,                  y+dt*k1,flowType,dataFlow,taup,mean_a_s,aa_s);
                k3 = FUNC_rhsSPARSER1D_inertial(time+3/4*dt,y+(1/4)*dt*k1+(1/4)*dt*k2,flowType,dataFlow,taup,mean_a_s,aa_s);
                y  = y+dt*((1/6)*k1+(1/6)*k2+(2/3)*k3);
            case 4 % RK4
                k1 = FUNC_rhsSPARSER1D_inertial(time       ,y           ,flowType,dataFlow,taup,mean_a_s,aa_s);
                k2 = FUNC_rhsSPARSER1D_inertial(time+0.5*dt,y+0.5*dt.*k1,flowType,dataFlow,taup,mean_a_s,aa_s);
                k3 = FUNC_rhsSPARSER1D_inertial(time+0.5*dt,y+0.5*dt.*k2,flowType,dataFlow,taup,mean_a_s,aa_s);
                k4 = FUNC_rhsSPARSER1D_inertial(time+dt    ,y+    dt.*k3,flowType,dataFlow,taup,mean_a_s,aa_s);
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
mean_xp = yy(0*nmp+1: 1*nmp,:);
mean_up = yy(1*nmp+1: 2*nmp,:);
xpxp    = yy(2*nmp+1: 3*nmp,:);
xpup    = yy(3*nmp+1: 4*nmp,:);
upup    = yy(4*nmp+1: 5*nmp,:);
axp     = yy(5*nmp+1: 6*nmp,:);
aup     = yy(6*nmp+1: 7*nmp,:);
S1.t    = t;

% Joining macro-particles -------------------------------------------------
w   = 1/nmp; % the weights would be all equal
for k=1:nt(1)
    smxp  = squeeze(mean_xp(:,k));
    smup  = squeeze(mean_up(:,k));
    sxpxp = squeeze(xpxp(:,k));
    sxpup = squeeze(xpup(:,k));
    supup = squeeze(upup(:,k));
    saxp  = squeeze(axp(:,k));
    saup  = squeeze(aup(:,k));
    S1.mean_xp(k,1) = sum(w*smxp);
    S1.mean_up(k,1) = sum(w*smup);
    S1.xpxp(k,1)    = sum(w*sxpxp) + sum(w*(smxp - S1.mean_xp(k)).^2);
    S1.xpup(k,1)    = sum(w*sxpup) + sum(w*(smxp - S1.mean_xp(k)).*(smup - S1.mean_up(k)));
    S1.upup(k,1)    = sum(w*supup) + sum(w*(smup - S1.mean_up(k)).^2);
    S1.axp(k,1)     = sum(w*saxp)  + sum(w*(mean_a_s - mean_a).*(smxp - S1.mean_xp(k)));
    S1.aup(k,1)     = sum(w*saup)  + sum(w*(mean_a_s - mean_a).*(smup - S1.mean_up(k)));
end

% % Post --------------------------------------------------------------------
% for j=1:nt(1)
%     for ix=1:nmpx
%         for iy=1:nmpy
%             Kxp                   = [S1.xpxp(ix,iy,j) S1.xpyp(ix,iy,j) ; S1.xpyp(ix,iy,j) S1.ypyp(ix,iy,j)];
%             S1.maxEigKxp(ix,iy,j) = max(eig(Kxp));
%             S1.FTLE(ix,iy,j)      = (1/(S1.t(j)-tLim(1)))*log(sqrt(S1.maxEigKxp(ix,iy,j)/S1.maxEigKxp(ix,iy,1)));
%             Kup                   = [S1.upup(ix,iy,j) S1.upvp(ix,iy,j) ; S1.upvp(ix,iy,j) S1.vpvp(ix,iy,j)];
%             S1.maxEigKup(ix,iy,j) = max(eig(Kup));
%             %             S1.FTLEv(ix,iy,j)     = (1/(S1.t(j)-tLim(1)))*log(sqrt(S1.maxEigKup(ix,iy,j)/S1.maxEigKup(ix,iy,1)));
%             S1.FTLEv(ix,iy,j)     = (1/(S1.t(j)-tLim(1)))*log(sqrt(S1.maxEigKup(ix,iy,j)));
%             S1.mu1(ix,iy,j)       = sqrt(S1.mean_xp(ix,iy,j).^2+S1.mean_yp(ix,iy,j).^2+S1.mean_up(ix,iy,j).^2+S1.mean_vp(ix,iy,j).^2);
%             K                     = [ ...
%                 S1.xpxp(ix,iy,j) S1.xpyp(ix,iy,j) S1.xpup(ix,iy,j) S1.xpvp(ix,iy,j) ;
%                 S1.xpyp(ix,iy,j) S1.ypyp(ix,iy,j) S1.ypup(ix,iy,j) S1.ypvp(ix,iy,j) ;
%                 S1.xpup(ix,iy,j) S1.ypup(ix,iy,j) S1.upup(ix,iy,j) S1.upvp(ix,iy,j) ;
%                 S1.xpvp(ix,iy,j) S1.ypvp(ix,iy,j) S1.upvp(ix,iy,j) S1.vpvp(ix,iy,j)];
%             S1.mu2(ix,iy,j)       = det(K);
%         end
%     end
% end
% S1.FTLE(:,:,1)  = 0;
% S1.FTLEv(:,:,1) = 0;
% S1.FTLE         = real(S1.FTLE);
% S1.FTLEv        = real(S1.FTLEv);

end
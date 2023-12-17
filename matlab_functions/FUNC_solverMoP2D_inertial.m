%% ========================================================================

function P1 = FUNC_solverMoP2D_inertial(xp0,yp0,up0,vp0,tLim,nt, ...
    flowType,dataFlow,f1Type,dataf1,timeMethod,taup)

%%

disp('FUNC_solverMoP_inertial')

% Small calculations ------------------------------------------------------
nt_2save = nt(1);
nt_2iter = nt(2);
Nt       = (nt(1)-1)*(nt(2)-1)+1;
dt       = (tLim(2)-tLim(1))/(Nt-1);
t        = linspace(tLim(1),tLim(2),nt_2save)';
ns       = length(dataf1.a1);

% Initializating ----------------------------------------------------------
npx  = length(xp0(:,1,1)); % clouds along x
npy  = length(xp0(1,:,1)); % clouds along y
np   = npx*npy; % total clouds
xp0c = xp0(:); % grid of initial particles put into columns
yp0c = yp0(:);
up0c = up0(:);
vp0c = vp0(:);
time = tLim(1);
yy   = [ xp0c ; yp0c ; up0c ; vp0c ];
y    = yy(:,1);

% ORDER ZERO --------------------------------------------------------------
for j=1:nt_2save-1            % loop used to save in time
    for dummyVar=1:nt_2iter-1 % loop that iterates in every time step to save
        switch timeMethod
            case 1 % EE
                k1 = FUNC_rhsMoP2D_o0_inertial(time       ,                        y,flowType,dataFlow,f1Type,dataf1,npx,npy,taup);
                y  = y+dt*k1;
            case 2 % RK2
                disp('not coded')
            case 3 % RK3 TVD
                k1 = FUNC_rhsMoP2D_o0_inertial(time       ,                        y,flowType,dataFlow,f1Type,dataf1,npx,npy,taup);
                k2 = FUNC_rhsMoP2D_o0_inertial(time+dt    ,                  y+dt*k1,flowType,dataFlow,f1Type,dataf1,npx,npy,taup);
                k3 = FUNC_rhsMoP2D_o0_inertial(time+3/4*dt,y+(1/4)*dt*k1+(1/4)*dt*k2,flowType,dataFlow,f1Type,dataf1,npx,npy,taup);
                y  = y+dt*((1/6)*k1+(1/6)*k2+(2/3)*k3);
            case 4 % RK4
                k1 = FUNC_rhsMoP2D_o0_inertial(time       ,y           ,flowType,dataFlow,f1Type,dataf1,npx,npy,taup);
                k2 = FUNC_rhsMoP2D_o0_inertial(time+0.5*dt,y+0.5*dt.*k1,flowType,dataFlow,f1Type,dataf1,npx,npy,taup);
                k3 = FUNC_rhsMoP2D_o0_inertial(time+0.5*dt,y+0.5*dt.*k2,flowType,dataFlow,f1Type,dataf1,npx,npy,taup);
                k4 = FUNC_rhsMoP2D_o0_inertial(time+dt    ,y+    dt.*k3,flowType,dataFlow,f1Type,dataf1,npx,npy,taup);
                y  = y+dt*(k1+2*k2+2*k3+k4)/6;
            otherwise
                disp('Error in timeMethod')
        end
        time = time + dt;
    end
    yy(:,j+1) = y; % Saving after iterating without saving
    if j==round(0.1*nt(1))
        fprintf('o0: 10%%..  t=%1.4f \n',time);
    elseif j==round(0.25*nt(1))
        fprintf('o0: 25%%..  t=%1.4f \n',time);
    elseif j==round(0.5*nt(1))
        fprintf('o0: 50%%..  t=%1.4f \n',time);
    elseif j==round(0.75*nt(1))
        fprintf('o0: 75%%..  t=%1.4f \n',time);
    elseif j==round(nt(1)-1)
        fprintf('o0: 100%%.  t=%1.4f \n',time);
    end
end
P1.o0.xp = reshape(yy(0*np+1:1*np,:),npx,npy,nt(1));
P1.o0.yp = reshape(yy(1*np+1:2*np,:),npx,npy,nt(1));
P1.o0.up = reshape(yy(2*np+1:3*np,:),npx,npy,nt(1));
P1.o0.vp = reshape(yy(3*np+1:4*np,:),npx,npy,nt(1));
P1.t     = t;

% ORDER EPSILON -----------------------------------------------------------
clear y yy
xp0c = xp0(:).*0; % grid of initial particles put into columns
yp0c = yp0(:).*0;
up0c = up0(:).*0;
vp0c = vp0(:).*0;
time = tLim(1);
yy   = [ xp0c ; yp0c ; up0c ; vp0c ];
y    = yy(:,1);
for j=1:nt_2save-1            % loop used to save in time
    for dummyVar=1:nt_2iter-1 % loop that iterates in every time step to save
        switch timeMethod
            case 1 % EE
                k1 = FUNC_rhsMoP2D_o1_inertial(time       ,                        y,flowType,dataFlow,f1Type,dataf1,npx,npy,taup,P1.o0,P1.t);
                y  = y+dt*k1;
            case 2 % RK2
                disp('not coded')
            case 3 % RK3 TVD
                k1 = FUNC_rhsMoP2D_o1_inertial(time       ,                        y,flowType,dataFlow,f1Type,dataf1,npx,npy,taup,P1.o0,P1.t);
                k2 = FUNC_rhsMoP2D_o1_inertial(time+dt    ,                  y+dt*k1,flowType,dataFlow,f1Type,dataf1,npx,npy,taup,P1.o0,P1.t);
                k3 = FUNC_rhsMoP2D_o1_inertial(time+3/4*dt,y+(1/4)*dt*k1+(1/4)*dt*k2,flowType,dataFlow,f1Type,dataf1,npx,npy,taup,P1.o0,P1.t);
                y  = y+dt*((1/6)*k1+(1/6)*k2+(2/3)*k3);
            case 4 % RK4
                k1 = FUNC_rhsMoP2D_o1_inertial(time       ,y           ,flowType,dataFlow,f1Type,dataf1,npx,npy,taup,P1.o0,P1.t);
                k2 = FUNC_rhsMoP2D_o1_inertial(time+0.5*dt,y+0.5*dt.*k1,flowType,dataFlow,f1Type,dataf1,npx,npy,taup,P1.o0,P1.t);
                k3 = FUNC_rhsMoP2D_o1_inertial(time+0.5*dt,y+0.5*dt.*k2,flowType,dataFlow,f1Type,dataf1,npx,npy,taup,P1.o0,P1.t);
                k4 = FUNC_rhsMoP2D_o1_inertial(time+dt    ,y+    dt.*k3,flowType,dataFlow,f1Type,dataf1,npx,npy,taup,P1.o0,P1.t);
                y  = y+dt*(k1+2*k2+2*k3+k4)/6;
            otherwise
                disp('Error in timeMethod')
        end
        time = time + dt;
    end
    yy(:,j+1) = y; % Saving after iterating without saving
    if j==round(0.1*nt(1))
        fprintf('o1: 10%%..  t=%1.4f \n',time);
    elseif j==round(0.25*nt(1))
        fprintf('o1: 25%%..  t=%1.4f \n',time);
    elseif j==round(0.5*nt(1))
        fprintf('o1: 50%%..  t=%1.4f \n',time);
    elseif j==round(0.75*nt(1))
        fprintf('o1: 75%%..  t=%1.4f \n',time);
    elseif j==round(nt(1)-1)
        fprintf('o1: 100%%.  t=%1.4f \n',time);
    end
end
P1.o1.xp = reshape(yy(0*np+1:1*np,:),npx,npy,nt(1));
P1.o1.yp = reshape(yy(1*np+1:2*np,:),npx,npy,nt(1));
P1.o1.up = reshape(yy(2*np+1:3*np,:),npx,npy,nt(1));
P1.o1.vp = reshape(yy(3*np+1:4*np,:),npx,npy,nt(1));

% ORDER EPSILON^2 ---------------------------------------------------------
clear y yy
xp0c = xp0(:).*0; % grid of initial particles put into columns
yp0c = yp0(:).*0;
up0c = up0(:).*0;
vp0c = vp0(:).*0;
time = tLim(1);
yy   = [ xp0c ; yp0c ; up0c ; vp0c ];
y    = yy(:,1);
for j=1:nt_2save-1            % loop used to save in time
    for dummyVar=1:nt_2iter-1 % loop that iterates in every time step to save
        switch timeMethod
            case 1 % EE
                k1 = FUNC_rhsMoP2D_o2_inertial(time       ,                        y,flowType,dataFlow,f1Type,dataf1,npx,npy,taup,P1.o0,P1.o1,P1.t);
                y  = y+dt*k1;
            case 2 % RK2
                disp('not coded')
            case 3 % RK3 TVD
                k1 = FUNC_rhsMoP2D_o2_inertial(time       ,                        y,flowType,dataFlow,f1Type,dataf1,npx,npy,taup,P1.o0,P1.o1,P1.t);
                k2 = FUNC_rhsMoP2D_o2_inertial(time+dt    ,                  y+dt*k1,flowType,dataFlow,f1Type,dataf1,npx,npy,taup,P1.o0,P1.o1,P1.t);
                k3 = FUNC_rhsMoP2D_o2_inertial(time+3/4*dt,y+(1/4)*dt*k1+(1/4)*dt*k2,flowType,dataFlow,f1Type,dataf1,npx,npy,taup,P1.o0,P1.o1,P1.t);
                y  = y+dt*((1/6)*k1+(1/6)*k2+(2/3)*k3);
            case 4 % RK4
                k1 = FUNC_rhsMoP2D_o2_inertial(time       ,y           ,flowType,dataFlow,f1Type,dataf1,npx,npy,taup,P1.o0,P1.o1,P1.t);
                k2 = FUNC_rhsMoP2D_o2_inertial(time+0.5*dt,y+0.5*dt.*k1,flowType,dataFlow,f1Type,dataf1,npx,npy,taup,P1.o0,P1.o1,P1.t);
                k3 = FUNC_rhsMoP2D_o2_inertial(time+0.5*dt,y+0.5*dt.*k2,flowType,dataFlow,f1Type,dataf1,npx,npy,taup,P1.o0,P1.o1,P1.t);
                k4 = FUNC_rhsMoP2D_o2_inertial(time+dt    ,y+    dt.*k3,flowType,dataFlow,f1Type,dataf1,npx,npy,taup,P1.o0,P1.o1,P1.t);
                y  = y+dt*(k1+2*k2+2*k3+k4)/6;
            otherwise
                disp('Error in timeMethod')
        end
        time = time + dt;
    end
    yy(:,j+1) = y; % Saving after iterating without saving
    if j==round(0.1*nt(1))
        fprintf('o2: 10%%..  t=%1.4f \n',time);
    elseif j==round(0.25*nt(1))
        fprintf('o2: 25%%..  t=%1.4f \n',time);
    elseif j==round(0.5*nt(1))
        fprintf('o2: 50%%..  t=%1.4f \n',time);
    elseif j==round(0.75*nt(1))
        fprintf('o2: 75%%..  t=%1.4f \n',time);
    elseif j==round(nt(1)-1)
        fprintf('o2: 100%%.  t=%1.4f \n',time);
    end
end
P1.o2.xp = reshape(yy(0*np+1:1*np,:),npx,npy,nt(1));
P1.o2.yp = reshape(yy(1*np+1:2*np,:),npx,npy,nt(1));
P1.o2.up = reshape(yy(2*np+1:3*np,:),npx,npy,nt(1));
P1.o2.vp = reshape(yy(3*np+1:4*np,:),npx,npy,nt(1));

% MONTE CARLO FAST SAMPLING -----------------------------------------------
prime_a1 = dataf1.a1-1;
for j=1:nt(1)
    for i=1:ns
        P1.xp(:,:,i,j) = P1.o0.xp(:,:,j) + prime_a1(i).*P1.o1.xp(:,:,j) + prime_a1(i).^2.*P1.o2.xp(:,:,j);
        P1.yp(:,:,i,j) = P1.o0.yp(:,:,j) + prime_a1(i).*P1.o1.yp(:,:,j) + prime_a1(i).^2.*P1.o2.yp(:,:,j);
        P1.up(:,:,i,j) = P1.o0.up(:,:,j) + prime_a1(i).*P1.o1.up(:,:,j) + prime_a1(i).^2.*P1.o2.up(:,:,j);
        P1.vp(:,:,i,j) = P1.o0.vp(:,:,j) + prime_a1(i).*P1.o1.vp(:,:,j) + prime_a1(i).^2.*P1.o2.vp(:,:,j);
    end
end

% Post MC -----------------------------------------------------------------
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
            %             T                     = abs(P1.t(j)-tLim(1));
            %             P1.FTLE(ix,iy,j)      = 1/T*log(sqrt(P1.maxEigKxp(ix,iy,j)/P1.maxEigKxp(ix,iy,1)));
            % Velocities
            Kup                   = [P1.upup(ix,iy,j) P1.upvp(ix,iy,j) ; P1.upvp(ix,iy,j) P1.vpvp(ix,iy,j)];
            if sum(isnan(Kup(:))) > 0 || sum(isinf(Kup(:))) > 0
                P1.maxEigKup(ix,iy,j) = 0;
            else
                P1.maxEigKup(ix,iy,j) = max(eig(Kxp));
            end
            % %             P1.FTLEv(ix,iy,j)     = (1/(P1.t(j)-tLim(1)))*log(sqrt(P1.maxEigKup(ix,iy,j)/P1.maxEigKup(ix,iy,1)));
            %             P1.FTLEv(ix,iy,j)     = 1/T*log(sqrt(P1.maxEigKup(ix,iy,j)));
            %             % Erros
            %             P1.mu1(ix,iy,j)       = sqrt(P1.mean_xp(ix,iy,j).^2+P1.mean_yp(ix,iy,j).^2+P1.mean_up(ix,iy,j).^2+P1.mean_vp(ix,iy,j).^2);
            %             K                     = [ ...
            %                 P1.xpxp(ix,iy,j) P1.xpyp(ix,iy,j) P1.xpup(ix,iy,j) P1.xpvp(ix,iy,j) ;
            %                 P1.xpyp(ix,iy,j) P1.ypyp(ix,iy,j) P1.ypup(ix,iy,j) P1.ypvp(ix,iy,j) ;
            %                 P1.xpup(ix,iy,j) P1.ypup(ix,iy,j) P1.upup(ix,iy,j) P1.upvp(ix,iy,j) ;
            %                 P1.xpvp(ix,iy,j) P1.ypvp(ix,iy,j) P1.upvp(ix,iy,j) P1.vpvp(ix,iy,j)];
            %             P1.mu2(ix,iy,j)       = det(K);
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
            %             P1.maxEigKa(ix,iy,j) = max(eig(Ka));
        end
    end
    %     j
end
% P1.FTLE(:,:,1) = 0;
% P1.FTLEv(:,:,1) = 0;

% for j=1:nt(1)
%     out.cx(:,j)  = [ out.xp(:,j) ; out.xp(1,j) ];
%     out.cy(:,j)  = [ out.yp(:,j) ; out.yp(1,j) ];
%     [out.curv_L(:,j),out.curv_R(:,j),out.curv_K(:,:,j)] = curvature([out.cx(:,j) out.cy(:,j)]);
%     [out.curv_maxR(j),out.curv_maxR_ind(j)]             = max(1./out.curv_R(:,j));
% end


% USING METHOD OF MOMENTS -------------------------------------------------
a1a1           = dataf1.a1a1;
a1a1a1         = dataf1.a1a1a1;
a1a1a1a1       = dataf1.a1a1a1a1;
P1.MoM.mean_xp = P1.o0.xp + a1a1*P1.o2.xp;
P1.MoM.mean_yp = P1.o0.yp + a1a1*P1.o2.yp;
P1.MoM.mean_up = P1.o0.up + a1a1*P1.o2.up;
P1.MoM.mean_vp = P1.o0.vp + a1a1*P1.o2.vp;
P1.MoM.xpxp    = a1a1*P1.o1.xp.^2        + a1a1a1*(P1.o2.xp.*P1.o1.xp + P1.o1.xp.*P1.o2.xp) + (3*a1a1^2+a1a1a1a1)*P1.o2.xp.^2;
P1.MoM.xpyp    = a1a1*P1.o1.xp.*P1.o1.yp + a1a1a1*(P1.o2.xp.*P1.o1.yp + P1.o1.xp.*P1.o2.yp) + (3*a1a1^2+a1a1a1a1)*P1.o2.xp.*P1.o2.yp;
P1.MoM.xpup    = a1a1*P1.o1.xp.*P1.o1.up + a1a1a1*(P1.o2.xp.*P1.o1.up + P1.o1.xp.*P1.o2.up) + (3*a1a1^2+a1a1a1a1)*P1.o2.xp.*P1.o2.up;
P1.MoM.xpvp    = a1a1*P1.o1.xp.*P1.o1.vp + a1a1a1*(P1.o2.xp.*P1.o1.vp + P1.o1.xp.*P1.o2.vp) + (3*a1a1^2+a1a1a1a1)*P1.o2.xp.*P1.o2.vp;
P1.MoM.ypyp    = a1a1*P1.o1.yp.^2        + a1a1a1*(P1.o2.yp.*P1.o1.yp + P1.o1.yp.*P1.o2.yp) + (3*a1a1^2+a1a1a1a1)*P1.o2.yp.^2;
P1.MoM.ypup    = a1a1*P1.o1.yp.*P1.o1.up + a1a1a1*(P1.o2.yp.*P1.o1.up + P1.o1.yp.*P1.o2.up) + (3*a1a1^2+a1a1a1a1)*P1.o2.yp.*P1.o2.up;
P1.MoM.ypvp    = a1a1*P1.o1.yp.*P1.o1.vp + a1a1a1*(P1.o2.yp.*P1.o1.vp + P1.o1.yp.*P1.o2.vp) + (3*a1a1^2+a1a1a1a1)*P1.o2.yp.*P1.o2.vp;
P1.MoM.upup    = a1a1*P1.o1.up.^2        + a1a1a1*(P1.o2.up.*P1.o1.up + P1.o1.up.*P1.o2.up) + (3*a1a1^2+a1a1a1a1)*P1.o2.up.^2;
P1.MoM.upvp    = a1a1*P1.o1.up.*P1.o1.vp + a1a1a1*(P1.o2.up.*P1.o1.vp + P1.o1.up.*P1.o2.vp) + (3*a1a1^2+a1a1a1a1)*P1.o2.up.*P1.o2.vp;
P1.MoM.vpvp    = a1a1*P1.o1.vp.^2        + a1a1a1*(P1.o2.vp.*P1.o1.vp + P1.o1.vp.*P1.o2.vp) + (3*a1a1^2+a1a1a1a1)*P1.o2.vp.^2;

P1.MoM.xpxpxp = - 3*a1a1^2*P1.o1.xp.^2.*P1.o2.xp - 3*a1a1^2*P1.o2.xp.^3 + 2*a1a1^3*P1.o2.xp.^3;
P1.MoM.upupup = - 3*a1a1^2*P1.o1.up.^2.*P1.o2.up - 3*a1a1^2*P1.o2.up.^3 + 2*a1a1^3*P1.o2.up.^3;

for j=1:nt(1)
    for ix=1:npx
        for iy=1:npy
            % Positions
            Kxp                       = [P1.MoM.xpxp(ix,iy,j) P1.MoM.xpyp(ix,iy,j) ; P1.MoM.xpyp(ix,iy,j) P1.MoM.ypyp(ix,iy,j)];
            if sum(isnan(Kxp(:))) > 0 || sum(isinf(Kxp(:))) > 0
                P1.MoM.maxEigKxp(ix,iy,j) = 0;
            else
                P1.MoM.maxEigKxp(ix,iy,j) = max(eig(Kxp));
            end
            % Velocities
            Kup                   = [P1.MoM.upup(ix,iy,j) P1.MoM.upvp(ix,iy,j) ; P1.MoM.upvp(ix,iy,j) P1.MoM.vpvp(ix,iy,j)];
            if sum(isnan(Kup(:))) > 0 || sum(isinf(Kup(:))) > 0
                P1.MoM.maxEigKup(ix,iy,j) = 0;
            else
                P1.MoM.maxEigKup(ix,iy,j) = max(eig(Kup));
            end
        end
    end
end




end









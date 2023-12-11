%% ========================================================================

function out = FUNC_rhsPSIC2D_inertial(time,y,flowType,dataFlow,f1Type,dataf1,npx,npy,taup)

% Variables ---------------------------------------------------------------
np = npx*npy;
xp = y(0*np+1:1*np);
yp = y(1*np+1:2*np);
up = y(2*np+1:3*np);
vp = y(3*np+1:4*np);

% Reading functions -------------------------------------------------------
% Flow
out1 = FUNC_flow(time,xp,yp,xp*0,flowType,dataFlow);
u    = out1.u;
v    = out1.v;

% Forcing function 
f1 = FUNC_f1_PSIC(u,v,u*0,up,vp,up*0,f1Type,dataf1);

% Equations ---------------------------------------------------------------
dxpdt = up;
dypdt = vp;
dupdt = (f1/taup).*(u - up);
dvpdt = (f1/taup).*(v - vp);

% Output ------------------------------------------------------------------
out = [ dxpdt ; dypdt ; dupdt ; dvpdt ];

end
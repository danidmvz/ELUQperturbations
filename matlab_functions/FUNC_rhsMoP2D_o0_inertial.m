%% ========================================================================

function out = FUNC_rhsMoP2D_o0_inertial(time,y,flowType,dataFlow,f1Type,dataf1,npx,npy,taup)
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

% Equations ---------------------------------------------------------------
dxpdt = up;
dypdt = vp;
dupdt = (dataf1.mean_a1/taup).*(u - up);
dvpdt = (dataf1.mean_a1/taup).*(v - vp);

% Output ------------------------------------------------------------------
out = [ dxpdt ; dypdt ; dupdt ; dvpdt ];

end
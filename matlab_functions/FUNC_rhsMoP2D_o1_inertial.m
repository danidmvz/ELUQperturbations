%% ========================================================================

function out = FUNC_rhsMoP2D_o1_inertial(time,y,flowType,dataFlow,f1Type,npx,npy,taup,o0,t)

% Variables ---------------------------------------------------------------
np = npx*npy;
nt = length(t);
xp = y(0*np+1:1*np);
yp = y(1*np+1:2*np);
up = y(2*np+1:3*np);
vp = y(3*np+1:4*np);

if flowType == 2
   mod = sqrt(xp.^2+yp.^2); 
   up  = up.* (mod>=1);
   vp  = vp.* (mod>=1);
end

% Order zero
xp00 = reshape(o0.xp,np,nt);
yp00 = reshape(o0.yp,np,nt);
up00 = reshape(o0.up,np,nt);
vp00 = reshape(o0.vp,np,nt);
% for i=1:np
%     xp0(i,1) = interp1(t,xp00(i,:)',time,'cubic');
%     yp0(i,1) = interp1(t,yp00(i,:)',time,'cubic');
%     up0(i,1) = interp1(t,up00(i,:)',time,'cubic');
%     vp0(i,1) = interp1(t,vp00(i,:)',time,'cubic');
% end
xp0 = interp1(t,xp00',time,'cubic'); xp0 = xp0';
yp0 = interp1(t,yp00',time,'cubic'); yp0 = yp0';
up0 = interp1(t,up00',time,'cubic'); up0 = up0';
vp0 = interp1(t,vp00',time,'cubic'); vp0 = vp0';

% Reading functions -------------------------------------------------------
% Flow
out1 = FUNC_flow(time,xp0,yp0,xp0*0,flowType,dataFlow);
u0    = out1.u;
v0    = out1.v;
dudx0 = out1.dudx;
dudy0 = out1.dudy;
dvdx0 = out1.dvdx;
dvdy0 = out1.dvdy;

u1 = xp.*dudx0 + yp.*dudy0;
v1 = xp.*dvdx0 + yp.*dvdy0;

% Forcing function
% f1 = FUNC_f1_PSIC(u,v,u*0,up,vp,up*0,f1Type,dataf1);

% Equations ---------------------------------------------------------------
dxpdt = up;
dypdt = vp;
dupdt = (1/taup).*( (u1-up) + (u0-up0) );
dvpdt = (1/taup).*( (v1-vp) + (v0-vp0) );

% Output ------------------------------------------------------------------
out = [ dxpdt ; dypdt ; dupdt ; dvpdt ];

end
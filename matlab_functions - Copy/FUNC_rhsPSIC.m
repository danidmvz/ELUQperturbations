%% RHS particle phase =====================================================

function out = FUNC_rhsPSIC(time,y,alpha1,alpha2,flowType,TfieldType,f1Type,f2Type,Re_inf,M_inf,St,Pr,cr,dp)

np = length(y)/7;

xp = y(0*np+1:1*np,1);
yp = y(1*np+1:2*np,1);
zp = y(2*np+1:3*np,1);
up = y(3*np+1:4*np,1);
vp = y(4*np+1:5*np,1);
wp = y(5*np+1:6*np,1);
Tp = y(6*np+1:7*np,1);

% Reading functions -------------------------------------------------------

% Flow
out1 = FUNC_flow(time,xp,yp,zp,flowType);
u = out1.u;
v = out1.v;
w = out1.w;

% Temperature field
if flowType ~= 100 % No Isotropic turbulence
    out1 = FUNC_Tfield(time,xp,yp,zp,TfieldType);
end
T = out1.T;

% Correction f1
if f1Type == 100 % Isotropic turbulence
    out1 = FUNC_f1_isoTurb(time,xp,yp,zp,u,v,w,up,vp,wp,Re_inf,M_inf,dp,f1Type);
else
    out1 = FUNC_f1(time,u,v,w,up,vp,wp,f1Type,Re_inf,dp);
end
f1 = alpha1.*out1.f1; clear out1

% Correction f2
if f2Type == 100 % Isotropic turbulence
    out1 = FUNC_f2_isoTurb(time,xp,yp,zp,u,v,w,up,vp,wp,Re_inf,Pr,dp,f2Type);
else
    out1 = FUNC_f2(time,u,v,w,up,vp,wp,f2Type,Re_inf,Pr,dp);
end
f2 = alpha2.*out1.f2;

% Equations ---------------------------------------------------------------

dxpdt = up;
dypdt = vp;
dzpdt = wp;

dupdt = (f1/St).*(u-up);
dvpdt = (f1/St).*(v-vp);
dwpdt = (f1/St).*(w-wp);

dTpdt = (2*cr/(3*Pr))*(f2/St).*(T-Tp);

% Output ------------------------------------------------------------------

out = [ ...
    dxpdt ; dypdt ; dzpdt ; ...
    dupdt ; dvpdt ; dwpdt ; ...
    dTpdt];

end
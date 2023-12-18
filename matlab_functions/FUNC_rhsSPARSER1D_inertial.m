%% ========================================================================

function out = FUNC_rhsSPARSER1D_inertial(time,y,flowType,dataFlow,taup,mean_a_s,aa_s)

% Variables ---------------------------------------------------------------
nmp     = length(mean_a_s);
mean_xp = y( 0*nmp+1: 1*nmp);
mean_up = y( 1*nmp+1: 2*nmp);
xpxp    = y( 2*nmp+1: 3*nmp);
xpup    = y( 3*nmp+1: 4*nmp);
upup    = y( 4*nmp+1: 5*nmp);
axp     = y( 5*nmp+1: 6*nmp);
aup     = y( 6*nmp+1: 7*nmp);

% Evaluation of the carrier flow derivatives ------------------------------
out1  = FUNC_flow(time,mean_xp,mean_xp*0,mean_xp*0,flowType,dataFlow);
u     = out1.u;
dudx  = out1.dudx;
d2udx = out1.d2udx;

% Closure -----------------------------------------------------------------
mean_u = u + 0.5*xpxp.*d2udx;
xpu    = xpxp.*dudx;
upu    = xpup.*dudx;
uu     =  xpu.*dudx;

% Definitions with the relative velocity ----------------------------------
xpax = xpu-xpup;
upax = upu-upup;
au   = axp.*dudx;
aax  = au-aup;

% Evaluation of the drag correction derivatives ---------------------------
% F = a_i*psi_i(a), but here F=a
f1        = mean_a_s;
df1dalpha = 1;
df1dax    = 0;
d2f1dax   = 0;
d2f1dadax = 0;
mean_f1   = f1 + 0.5*( aax.*d2f1dadax + (uu -2*upu + upup).*d2f1dax );
f1a       = aa_s.*df1dalpha;

% Closure terms involving f1  ---------------------------------------------
xpf1 = axp.*df1dalpha + xpax.*df1dax;
upf1 = aup.*df1dalpha + upax.*df1dax;
f1u  = au;

% EQUATIONS ===============================================================
ddt_mean_xp = mean_up;
ddt_mean_up = (1/taup)*( mean_f1.*(mean_u-mean_up) + f1u - upf1 );
ddt_xpxp    = 2*xpup;
ddt_xpup    = upup + ( mean_f1.*(xpu-xpup) + xpf1.*(mean_u-mean_up) )/taup;
ddt_upup    = ( mean_f1.*(2*upu-2*upup) + 2*upf1.*(mean_u-mean_up) )*(1/taup);
ddt_axp     = aup;
ddt_aup     = ( mean_f1.*(au-aup) +f1a.*(mean_u-mean_up) )/taup;

% Outputs -----------------------------------------------------------------
out = [ddt_mean_xp ; ddt_mean_up ; ddt_xpxp ; ddt_xpup ; ddt_upup ; ...
    ddt_axp ; ddt_aup];

% Checks ------------------------------------------------------------------
% if sum(imag(xpxp)) > 0 || sum(imag(xpyp)) > 0 || sum(imag(ypyp)) > 0
%     disp('SPARSE seocnd moments imag')
% end

% if sum(xpxp < 0) || sum(ypyp < 0)
%     disp('SPARSE xpxp or ypyp negative')
% end
    
end




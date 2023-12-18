%% ========================================================================


function out = FUNC_rhsSPARSER2D_inertial(time,y,flowType,dataFlow,taup,mean_a_s,aa_s,npx,npy,M)

% Variables ---------------------------------------------------------------
NP      = npx*npy*M;
mean_xp = y( 0*NP+1: 1*NP);
mean_yp = y( 1*NP+1: 2*NP);
mean_up = y( 2*NP+1: 3*NP);
mean_vp = y( 3*NP+1: 4*NP);
xpxp    = y( 4*NP+1: 5*NP);
xpyp    = y( 5*NP+1: 6*NP);
xpup    = y( 6*NP+1: 7*NP);
xpvp    = y( 7*NP+1: 8*NP);
ypyp    = y( 8*NP+1: 9*NP);
ypup    = y( 9*NP+1:10*NP);
ypvp    = y(10*NP+1:11*NP);
upup    = y(11*NP+1:12*NP);
upvp    = y(12*NP+1:13*NP);
vpvp    = y(13*NP+1:14*NP);
axp     = y(14*NP+1:15*NP);
ayp     = y(15*NP+1:16*NP);
aup     = y(16*NP+1:17*NP);
avp     = y(17*NP+1:18*NP);
aa      = aa_s(:);
mean_a  = mean_a_s(:);

% Evaluation of the carrier flow derivatives ------------------------------
out1 = FUNC_flow(time,mean_xp,mean_yp,mean_xp*0,flowType,dataFlow);
% For u
u       = out1.u;
dudx    = out1.dudx;
dudy    = out1.dudy;
d2udx   = out1.d2udx;
d2udy   = out1.d2udy;
d2udxdy = out1.d2udxdy;
% For v
v       = out1.v;
dvdx    = out1.dvdx;
dvdy    = out1.dvdy;
d2vdx   = out1.d2vdx;
d2vdy   = out1.d2vdy;
d2vdxdy = out1.d2vdxdy;

% Closure -----------------------------------------------------------------
mean_u = u + 0.5*(xpxp.*d2udx + 2*xpyp.*d2udxdy + ypyp.*d2udy);
mean_v = v + 0.5*(xpxp.*d2vdx + 2*xpyp.*d2vdxdy + ypyp.*d2vdy);

% Terms xpu
xpu = xpxp.*dudx + xpyp.*dudy;
xpv = xpxp.*dvdx + xpyp.*dvdy;
ypu = xpyp.*dudx + ypyp.*dudy;
ypv = xpyp.*dvdx + ypyp.*dvdy;

% Terms upu
upu = xpup.*dudx + ypup.*dudy;
upv = xpup.*dvdx + ypup.*dvdy;
vpu = xpvp.*dudx + ypvp.*dudy;
vpv = xpvp.*dvdx + ypvp.*dvdy;

% Terms uu
uu = xpu.*dudx + ypu.*dudy;
vv = xpv.*dvdx + ypv.*dvdy;
uv = xpv.*dudx + ypv.*dudy;

% Evaluation of the drag correction derivatives ---------------------------
f1         = mean_a;
df1dalpha  = 1;
df1dax     = 0;
df1day     = 0;
d2f1dax    = 0;
d2f1day    = 0;
d2f1daxday = 0;

mean_f1 = f1 + 0.5*( ...
    + (uu -2*upu + upup).*d2f1dax ...
    + (vv -2*vpv + vpvp).*d2f1day ...
    + 2*(uv - vpu - upv + upvp).*d2f1daxday);

% Definitions with the relative velocity ----------------------------------
% Terms xpa
xpax = xpu-xpup;
xpay = xpv-xpvp;
ypax = ypu-ypup;
ypay = ypv-ypvp;
% Terms upa
upax = upu-upup;
upay = upv-upvp;
vpax = vpu-upvp;
vpay = vpv-vpvp;

% Closure terms involving f1 and f2 ---------------------------------------
xpf1 = axp.*df1dalpha + xpax.*df1dax + xpay.*df1day;
ypf1 = ayp.*df1dalpha + ypax.*df1dax + ypay.*df1day;
upf1 = aup.*df1dalpha + upax.*df1dax + upay.*df1day;
vpf1 = avp.*df1dalpha + vpax.*df1dax + vpay.*df1day;

f1u = xpf1.*dudx + ypf1.*dudy;
f1v = xpf1.*dvdx + ypf1.*dvdy;
au  = axp.*dudx  + ayp.*dudy;
av  = axp.*dvdx  + ayp.*dvdy;
f1a = aa.*df1dalpha;


% EQUATIONS ===============================================================
% Means
ddt_mean_xp = mean_up;
ddt_mean_yp = mean_vp;
ddt_mean_up = (1/taup)*( mean_f1.*(mean_u-mean_up) + f1u - upf1 );
ddt_mean_vp = (1/taup)*( mean_f1.*(mean_v-mean_vp) + f1v - vpf1 );

% Terms xpxp
ddt_xpxp = 2*xpup;
ddt_ypyp = 2*ypvp;
ddt_xpyp = xpvp+ypup;

% Terms xpup
ddt_xpup = upup + ( mean_f1.*(xpu-xpup) + xpf1.*(mean_u-mean_up) )/taup;
ddt_xpvp = upvp + ( mean_f1.*(xpv-xpvp) + xpf1.*(mean_v-mean_vp) )/taup;
ddt_ypup = upvp + ( mean_f1.*(ypu-ypup) + ypf1.*(mean_u-mean_up) )/taup;
ddt_ypvp = vpvp + ( mean_f1.*(ypv-ypvp) + ypf1.*(mean_v-mean_vp) )/taup;

% Terms upup
ddt_upup = ( mean_f1.*(2*upu-2*upup) + 2*upf1.*(mean_u-mean_up) )*(1/taup);
ddt_vpvp = ( mean_f1.*(2*vpv-2*vpvp) + 2*vpf1.*(mean_v-mean_vp) )*(1/taup);
ddt_upvp = ( mean_f1.*(vpu+upv-2*upvp) + vpf1.*(mean_u-mean_up) + upf1.*(mean_v-mean_vp) )/taup;

% Terms coeff with particle phase
ddt_axp = aup;
ddt_ayp = avp;
ddt_aup = (mean_f1.*(au-aup)+f1a.*(mean_u-mean_up))/taup;
ddt_avp = (mean_f1.*(av-avp)+f1a.*(mean_v-mean_vp))/taup;

% Outputs -----------------------------------------------------------------
out = [ ...
    ddt_mean_xp ; ddt_mean_yp ; ddt_mean_up ; ddt_mean_vp ; ...
    ddt_xpxp    ; ddt_xpyp    ; ddt_xpup    ; ddt_xpvp    ; ...
    ddt_ypyp    ; ddt_ypup    ; ddt_ypvp    ; ...
    ddt_upup    ; ddt_upvp    ; ...
    ddt_vpvp    ; ...
    ddt_axp     ; ddt_ayp     ; ddt_aup     ; ddt_avp];

   
end




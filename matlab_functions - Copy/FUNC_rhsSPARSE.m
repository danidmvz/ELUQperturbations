%% RHS MoM 3D ==========================================================

function out = FUNC_rhsSPARSE(time,y,mean_alpha1,mean_alpha2,cm2_alpha1,alpha1alpha2,flowType,TfieldType,f1Type,f2Type,Re_inf,M_inf,St,Pr,cr,dp,nmp)

% Variables ---------------------------------------------------------------
% Means
mean_xp = y( 0*nmp+1: 1*nmp);
mean_yp = y( 1*nmp+1: 2*nmp);
mean_zp = y( 2*nmp+1: 3*nmp);
mean_up = y( 3*nmp+1: 4*nmp);
mean_vp = y( 4*nmp+1: 5*nmp);
mean_wp = y( 5*nmp+1: 6*nmp);
mean_Tp = y( 6*nmp+1: 7*nmp);

% xp with xp
cm2_xp  = y( 7*nmp+1: 8*nmp);
cm2_yp  = y( 8*nmp+1: 9*nmp);
cm2_zp  = y( 9*nmp+1:10*nmp);
xpyp    = y(10*nmp+1:11*nmp);
xpzp    = y(11*nmp+1:12*nmp);
ypzp    = y(12*nmp+1:13*nmp);

% xp with up
xpup    = y(13*nmp+1:14*nmp);
xpvp    = y(14*nmp+1:15*nmp);
xpwp    = y(15*nmp+1:16*nmp);
ypup    = y(16*nmp+1:17*nmp);
ypvp    = y(17*nmp+1:18*nmp);
ypwp    = y(18*nmp+1:19*nmp);
zpup    = y(19*nmp+1:20*nmp);
zpvp    = y(20*nmp+1:21*nmp);
zpwp    = y(21*nmp+1:22*nmp);

% xp with Tp
xpTp    = y(22*nmp+1:23*nmp);
ypTp    = y(23*nmp+1:24*nmp);
zpTp    = y(24*nmp+1:25*nmp);

% up with up
cm2_up  = y(25*nmp+1:26*nmp);
cm2_vp  = y(26*nmp+1:27*nmp);
cm2_wp  = y(27*nmp+1:28*nmp);
upvp    = y(28*nmp+1:29*nmp);
upwp    = y(29*nmp+1:30*nmp);
vpwp    = y(30*nmp+1:31*nmp);

% up with Tp
upTp    = y(31*nmp+1:32*nmp);
vpTp    = y(32*nmp+1:33*nmp);
wpTp    = y(33*nmp+1:34*nmp);

% Tp with Tp
cm2_Tp  = y(34*nmp+1:35*nmp);

% alpha1 with xp
xpalpha1 = y(35*nmp+1:36*nmp);
ypalpha1 = y(36*nmp+1:37*nmp);
zpalpha1 = y(37*nmp+1:38*nmp);

% alpha1 with up
upalpha1 = y(38*nmp+1:39*nmp);
vpalpha1 = y(39*nmp+1:40*nmp);
wpalpha1 = y(40*nmp+1:41*nmp);

% alpha2 with xp
xpalpha2 = y(41*nmp+1:42*nmp);
ypalpha2 = y(42*nmp+1:43*nmp);
zpalpha2 = y(43*nmp+1:44*nmp);

% alpha2 with up
upalpha2 = y(44*nmp+1:45*nmp);
vpalpha2 = y(45*nmp+1:46*nmp);
wpalpha2 = y(46*nmp+1:47*nmp);

% alpha1 with Tp
Tpalpha1 = y(47*nmp+1:48*nmp);

% alpha2 with Tp
Tpalpha2 = y(48*nmp+1:49*nmp);


% Evaluation of the carrier flow derivatives ------------------------------
out1        = FUNC_flow(time,mean_xp,mean_yp,mean_zp,flowType);
% For u
u_at_xp     = out1.u;
dudx_at_xp  = out1.dudx;
d2udx_at_xp = out1.d2udx;
dudy_at_xp  = out1.dudy;
d2udy_at_xp = out1.d2udy;
dudz_at_xp  = out1.dudz;
d2udz_at_xp = out1.d2udz;
% For v
v_at_xp     = out1.v;
dvdx_at_xp  = out1.dvdx;
d2vdx_at_xp = out1.d2vdx;
dvdy_at_xp  = out1.dvdy;
d2vdy_at_xp = out1.d2vdy;
dvdz_at_xp  = out1.dvdz;
d2vdz_at_xp = out1.d2vdz;
% For w
w_at_xp     = out1.w;
dwdx_at_xp  = out1.dwdx;
d2wdx_at_xp = out1.d2wdx;
dwdy_at_xp  = out1.dwdy;
d2wdy_at_xp = out1.d2wdy;
dwdz_at_xp  = out1.dwdz;
d2wdz_at_xp = out1.d2wdz;


% Evaluation of the carrier T field derivatives ---------------------------
if flowType ~= 100 % No Isotropic turbulence
    out1 = FUNC_Tfield(time,mean_xp,mean_yp,mean_zp,TfieldType);
end
T_at_xp     = out1.T;
dTdx_at_xp  = out1.dTdx;
d2Tdx_at_xp = out1.d2Tdx;
dTdy_at_xp  = out1.dTdy;
d2Tdy_at_xp = out1.d2Tdy;
dTdz_at_xp  = out1.dTdz;
d2Tdz_at_xp = out1.d2Tdz;


% Closure -----------------------------------------------------------------
% Mean u and T
mean_u = u_at_xp + 0.5*cm2_xp.*d2udx_at_xp + 0.5*cm2_yp.*d2udy_at_xp + 0.5*cm2_zp.*d2udz_at_xp;
mean_v = v_at_xp + 0.5*cm2_xp.*d2vdx_at_xp + 0.5*cm2_yp.*d2vdy_at_xp + 0.5*cm2_zp.*d2vdz_at_xp;
mean_w = w_at_xp + 0.5*cm2_xp.*d2wdx_at_xp + 0.5*cm2_yp.*d2wdy_at_xp + 0.5*cm2_zp.*d2wdz_at_xp;
mean_T = T_at_xp + 0.5*cm2_xp.*d2Tdx_at_xp + 0.5*cm2_yp.*d2Tdy_at_xp + 0.5*cm2_zp.*d2Tdz_at_xp;

% Terms xpu
xpu = cm2_xp.*dudx_at_xp +   xpyp.*dudy_at_xp +   xpzp.*dudz_at_xp;
xpv = cm2_xp.*dvdx_at_xp +   xpyp.*dvdy_at_xp +   xpzp.*dvdz_at_xp;
xpw = cm2_xp.*dwdx_at_xp +   xpyp.*dwdy_at_xp +   xpzp.*dwdz_at_xp;
ypu =   xpyp.*dudx_at_xp + cm2_yp.*dudy_at_xp +   ypzp.*dudz_at_xp;
ypv =   xpyp.*dvdx_at_xp + cm2_yp.*dvdy_at_xp +   ypzp.*dvdz_at_xp;
ypw =   xpyp.*dwdx_at_xp + cm2_yp.*dwdy_at_xp +   ypzp.*dwdz_at_xp;
zpu =   xpzp.*dudx_at_xp +   ypzp.*dudy_at_xp + cm2_zp.*dudz_at_xp;
zpv =   xpzp.*dvdx_at_xp +   ypzp.*dvdy_at_xp + cm2_zp.*dvdz_at_xp;
zpw =   xpzp.*dwdx_at_xp +   ypzp.*dwdy_at_xp + cm2_zp.*dwdz_at_xp;

% Terms upu
upu = xpup.*dudx_at_xp + ypup.*dudy_at_xp + zpup.*dudz_at_xp;
upv = xpup.*dvdx_at_xp + ypup.*dvdy_at_xp + zpup.*dvdz_at_xp;
upw = xpup.*dwdx_at_xp + ypup.*dwdy_at_xp + zpup.*dwdz_at_xp;
vpu = xpvp.*dudx_at_xp + ypvp.*dudy_at_xp + zpvp.*dudz_at_xp;
vpv = xpvp.*dvdx_at_xp + ypvp.*dvdy_at_xp + zpvp.*dvdz_at_xp;
vpw = xpvp.*dwdx_at_xp + ypvp.*dwdy_at_xp + zpvp.*dwdz_at_xp;
wpu = xpwp.*dudx_at_xp + ypwp.*dudy_at_xp + zpwp.*dudz_at_xp;
wpv = xpwp.*dvdx_at_xp + ypwp.*dvdy_at_xp + zpwp.*dvdz_at_xp;
wpw = xpwp.*dwdx_at_xp + ypwp.*dwdy_at_xp + zpwp.*dwdz_at_xp;

% % Terms uu
% cm2_u = xpu.*dudx_at_xp + ypu.*dudy_at_xp + ypw.*dvdz_at_xp;
% cm2_v = xpv.*dvdx_at_xp + ypv.*dvdy_at_xp + ypw.*dvdz_at_xp;
% cm2_w = xpw.*dwdx_at_xp + ypw.*dwdy_at_xp + zpw.*dwdz_at_xp;
%
% uv    = xpv.*dudx_at_xp + ypv.*dudy_at_xp + zpv.*dudz_at_xp;
% uw    = xpw.*dudx_at_xp + ypw.*dudy_at_xp + zpw.*dudz_at_xp;
% vw    = xpw.*dvdx_at_xp + ypw.*dvdy_at_xp + zpw.*dvdz_at_xp;

% Terms xpT
xpT   = cm2_xp.*dTdx_at_xp +   xpyp.*dTdy_at_xp +   xpzp.*dTdz_at_xp;
ypT   =   xpyp.*dTdx_at_xp + cm2_yp.*dTdy_at_xp +   ypzp.*dTdz_at_xp;
zpT   =   xpzp.*dTdx_at_xp +   ypzp.*dTdy_at_xp + cm2_zp.*dTdz_at_xp;

% Terms upT
upT   = xpup.*dTdx_at_xp + ypup.*dTdy_at_xp + zpup.*dTdz_at_xp;
vpT   = xpvp.*dTdx_at_xp + ypvp.*dTdy_at_xp + zpvp.*dTdz_at_xp;
wpT   = xpwp.*dTdx_at_xp + ypwp.*dTdy_at_xp + zpwp.*dTdz_at_xp;

% % Terms uT
% uT   = xpT.*dudx_at_xp + ypT.*dudy_at_xp + zpT.*dudz_at_xp;
% vT   = xpT.*dvdx_at_xp + ypT.*dvdy_at_xp + zpT.*dvdz_at_xp;
% wT   = xpT.*dwdx_at_xp + ypT.*dwdy_at_xp + zpT.*dwdz_at_xp;

% Term cm2_T (not needed)
% cm2_T = xpT.*dTdx_at_xp + ypT.*dTdy_at_xp + zpT.*dTdz_at_xp;

% Terms Tpu
Tpu = xpTp.*dudx_at_xp + ypTp.*dudy_at_xp + zpTp.*dudz_at_xp;
Tpv = xpTp.*dvdx_at_xp + ypTp.*dvdy_at_xp + zpTp.*dvdz_at_xp;
Tpw = xpTp.*dwdx_at_xp + ypTp.*dwdy_at_xp + zpTp.*dwdz_at_xp;

% Term TpT
TpT = xpTp.*dTdx_at_xp + ypTp.*dTdy_at_xp + zpTp.*dTdz_at_xp;

% Evaluation of the drag correction derivatives ---------------------------
if f1Type == 100 % Isotropic turbulence
    out1 = FUNC_f1_isoTurb(time,mean_xp,mean_yp,mean_zp,mean_u,mean_v,mean_w,mean_up,mean_vp,mean_wp,Re_inf,M_inf,dp,f1Type);
else
    out1 = FUNC_f1(time,mean_u,mean_v,mean_w,mean_up,mean_vp,mean_wp,f1Type,Re_inf,dp);
end
mean_f1      = mean_alpha1.*out1.f1;
mean_df1dax  = mean_alpha1.*out1.df1dax;
mean_df1day  = mean_alpha1.*out1.df1day;
mean_df1daz  = mean_alpha1.*out1.df1daz;
% mean_d2f1dax = out1.d2f1dax;
% mean_d2f1day = out1.d2f1day;
% mean_d2f1daz = out1.d2f1daz;


% Evaluation of the Nu correction derivatives -----------------------------
if f1Type == 100 % Isotropic turbulence
    out1 = FUNC_f2_isoTurb(time,mean_xp,mean_yp,mean_zp,mean_u,mean_v,mean_w,mean_up,mean_vp,mean_wp,Re_inf,Pr,dp,f2Type);
else
    out1 = FUNC_f2(time,mean_u,mean_v,mean_w,mean_up,mean_vp,mean_wp,f2Type,Re_inf,Pr,dp);
end
mean_f2      = mean_alpha2.*out1.f2;
mean_df2dax  = mean_alpha2.*out1.df2dax;
mean_df2day  = mean_alpha2.*out1.df2day;
mean_df2daz  = mean_alpha2.*out1.df2daz;
% mean_d2f2dax = out1.d2f2dax;
% mean_d2f2day = out1.d2f2day;
% mean_d2f2daz = out1.d2f2daz;

% % Definitions with the relative velocity ----------------------------------
% % Mean a
% mean_ax = mean_u-mean_up;
% mean_ay = mean_v-mean_vp;
% mean_az = mean_w-mean_wp;

% % Terms aa
% cm2_ax = cm2_u+cm2_up-2*upu;
% cm2_ay = cm2_v+cm2_vp-2*vpv;
% cm2_az = cm2_w+cm2_wp-2*wpw;
% axay   = uv+upvp-upv-vpu;
% axaz   = uw+upwp-upw-wpu;
% ayaz   = vw+vpwp-ypw-wpv;

% Terms xpa
xpax = xpu-xpup;
xpay = xpv-xpvp;
xpaz = xpw-xpwp;
ypax = ypu-ypup;
ypay = ypv-ypvp;
ypaz = ypw-ypwp;
zpax = zpu-zpup;
zpay = zpv-zpvp;
zpaz = zpw-zpwp;

% Terms upa
upax = upu-cm2_up;
upay = upv-upvp;
upaz = upw-upwp;
vpax = vpu-upvp;
vpay = vpv-cm2_vp;
vpaz = vpw-vpwp;
wpax = wpu-upwp;
wpay = wpv-vpwp;
wpaz = wpw-cm2_wp;

% Terms Tpa
Tpax = Tpu-upTp;
Tpay = Tpv-vpTp;
Tpaz = Tpw-wpTp;


% % Definitions with the relative temperature -------------------------------
% mean_b = mean_T-mean_Tp;
%
% % Terms ab
% axb = uT+upTp-upT-Tpu;
% ayb = vT+vpTp-vpT-Tpv;
% azb = wT+wpTp-wpT-Tpw;
%
% % Terms xpb
% xpb = xpT-xpTp;
% ypb = ypT-ypTp;
% zpb = zpT-zpTp;
%
% % Terms upb
% upb = upT-upTp;
% vpb = vpT-vpTp;
% wpb = wpT-wpTp;
%
% % Terms Tpb
% Tpb = TpT-cm2_Tp;

% Closure terms involving f1 and f2 ---------------------------------------
df1dalpha1 = mean_f1;
df2dalpha2 = mean_f2;
xpf1 = xpalpha1.*df1dalpha1+xpax.*mean_df1dax+xpay.*mean_df1day+xpaz.*mean_df1daz;
ypf1 = ypalpha1.*df1dalpha1+ypax.*mean_df1dax+ypay.*mean_df1day+ypaz.*mean_df1daz;
zpf1 = zpalpha1.*df1dalpha1+zpax.*mean_df1dax+zpay.*mean_df1day+zpaz.*mean_df1daz;
upf1 = upalpha1.*df1dalpha1+upax.*mean_df1dax+upay.*mean_df1day+upaz.*mean_df1daz;
vpf1 = vpalpha1.*df1dalpha1+vpax.*mean_df1dax+vpay.*mean_df1day+vpaz.*mean_df1daz;
wpf1 = wpalpha1.*df1dalpha1+wpax.*mean_df1dax+wpay.*mean_df1day+wpaz.*mean_df1daz;

xpf2 = xpalpha2.*df2dalpha2+xpax.*mean_df2dax+xpay.*mean_df2day+xpaz.*mean_df2daz;
ypf2 = ypalpha2.*df2dalpha2+ypax.*mean_df2dax+ypay.*mean_df2day+ypaz.*mean_df2daz;
zpf2 = zpalpha2.*df2dalpha2+zpax.*mean_df2dax+zpay.*mean_df2day+zpaz.*mean_df2daz;
upf2 = upalpha2.*df2dalpha2+upax.*mean_df2dax+upay.*mean_df2day+upaz.*mean_df2daz;
vpf2 = vpalpha2.*df2dalpha2+vpax.*mean_df2dax+vpay.*mean_df2day+vpaz.*mean_df2daz;
wpf2 = wpalpha2.*df2dalpha2+wpax.*mean_df2dax+wpay.*mean_df2day+wpaz.*mean_df2daz;


Tpf1 = Tpalpha1.*df1dalpha1+Tpax.*mean_df1dax+Tpay.*mean_df1day+Tpaz.*mean_df1daz;
Tpf2 = Tpalpha2.*df2dalpha2+Tpax.*mean_df2dax+Tpay.*mean_df2day+Tpaz.*mean_df2daz;

f1u = xpf1.*dudx_at_xp + ypf1.*dudy_at_xp + zpf1.*dudz_at_xp;
f1v = xpf1.*dvdx_at_xp + ypf1.*dvdy_at_xp + zpf1.*dvdz_at_xp;
f1w = xpf1.*dwdx_at_xp + ypf1.*dwdy_at_xp + zpf1.*dwdz_at_xp;

f2T = xpf2.*dTdx_at_xp + ypf2.*dTdy_at_xp + zpf2.*dTdz_at_xp;

ualpha1 = xpalpha1.*dudx_at_xp + ypalpha1.*dudy_at_xp + zpalpha1.*dudz_at_xp;
valpha1 = xpalpha1.*dvdx_at_xp + ypalpha1.*dvdy_at_xp + zpalpha1.*dvdz_at_xp;
walpha1 = xpalpha1.*dwdx_at_xp + ypalpha1.*dwdy_at_xp + zpalpha1.*dwdz_at_xp;
Talpha1 = xpalpha1.*dTdx_at_xp + ypalpha1.*dTdy_at_xp + zpalpha1.*dTdz_at_xp;

ualpha2 = xpalpha2.*dudx_at_xp + ypalpha2.*dudy_at_xp + zpalpha2.*dudz_at_xp;
valpha2 = xpalpha2.*dvdx_at_xp + ypalpha2.*dvdy_at_xp + zpalpha2.*dvdz_at_xp;
walpha2 = xpalpha2.*dwdx_at_xp + ypalpha2.*dwdy_at_xp + zpalpha2.*dwdz_at_xp;
Talpha2 = xpalpha2.*dTdx_at_xp + ypalpha2.*dTdy_at_xp + zpalpha2.*dTdz_at_xp;

alpha1ax = ualpha1-upalpha1;
alpha1ay = valpha1-vpalpha1;
alpha1az = walpha1-wpalpha1;
alpha2ax = ualpha2-upalpha2;
alpha2ay = valpha2-vpalpha2;
alpha2az = walpha2-wpalpha2;
f1alpha1 =   cm2_alpha1.*df1dalpha1 + alpha1ax.*mean_df1dax+alpha1ay.*mean_df1day+alpha1az.*mean_df1daz;
f1alpha2 = alpha1alpha2.*df1dalpha1 + alpha2ax.*mean_df1dax+alpha2ay.*mean_df1day+alpha2az.*mean_df1daz;

% EQUATIONS ===============================================================

% Means
ddt_mean_xp = mean_up;
ddt_mean_yp = mean_vp;
ddt_mean_zp = mean_wp;

ddt_mean_up = ( mean_f1.*(mean_u-mean_up) + f1u - upf1 )/St;
ddt_mean_vp = ( mean_f1.*(mean_v-mean_vp) + f1v - vpf1 )/St;
ddt_mean_wp = ( mean_f1.*(mean_w-mean_wp) + f1w - wpf1 )/St;

ddt_mean_Tp = ( mean_f2.*(mean_T-mean_Tp) + f2T - Tpf2 )/St;

% Terms xpxp
ddt_cm2_xp = 2*xpup;
ddt_cm2_yp = 2*ypvp;
ddt_cm2_zp = 2*zpwp;

ddt_xpyp   = xpvp+ypup;
ddt_xpzp   = xpwp+zpup;
ddt_ypzp   = ypwp+zpvp;

% Terms xpup
ddt_xpup = cm2_up + ( mean_f1.*(xpu-xpup) + xpf1.*(mean_u-mean_up) )/St;
ddt_xpvp =   upvp + ( mean_f1.*(xpv-xpvp) + xpf1.*(mean_v-mean_vp) )/St;
ddt_xpwp =   upwp + ( mean_f1.*(xpw-xpwp) + xpf1.*(mean_w-mean_wp) )/St;

ddt_ypup =   upvp + ( mean_f1.*(ypu-ypup) + ypf1.*(mean_u-mean_up) )/St;
ddt_ypvp = cm2_vp + ( mean_f1.*(ypv-ypvp) + ypf1.*(mean_v-mean_vp) )/St;
ddt_ypwp =   vpwp + ( mean_f1.*(ypw-ypwp) + ypf1.*(mean_w-mean_wp) )/St;

ddt_zpup =   upwp + ( mean_f1.*(zpu-zpup) + zpf1.*(mean_u-mean_up) )/St;
ddt_zpvp =   vpwp + ( mean_f1.*(zpv-zpvp) + zpf1.*(mean_v-mean_vp) )/St;
ddt_zpwp = cm2_wp + ( mean_f1.*(zpw-zpwp) + zpf1.*(mean_w-mean_wp) )/St;

% Terms xpTp
ddt_xpTp = upTp+(2*cr/(3*Pr*St))*( mean_f2.*(xpT-xpTp) + xpf2.*(mean_T-mean_Tp) );
ddt_ypTp = vpTp+(2*cr/(3*Pr*St))*( mean_f2.*(ypT-ypTp) + ypf2.*(mean_T-mean_Tp) );
ddt_zpTp = wpTp+(2*cr/(3*Pr*St))*( mean_f2.*(zpT-zpTp) + zpf2.*(mean_T-mean_Tp) );

% Terms upup
ddt_cm2_up = ( mean_f1.*(2*upu-2*cm2_up) + 2*upf1.*(mean_u-mean_up) )*(1/St);
ddt_cm2_vp = ( mean_f1.*(2*vpv-2*cm2_vp) + 2*vpf1.*(mean_v-mean_vp) )*(1/St);
ddt_cm2_wp = ( mean_f1.*(2*wpw-2*cm2_wp) + 2*wpf1.*(mean_w-mean_wp) )*(1/St);

ddt_upvp   = ( mean_f1.*(vpu+upv-2*upvp) + vpf1.*(mean_u-mean_up) + upf1.*(mean_v-mean_vp) )/St;
ddt_upwp   = ( mean_f1.*(wpu+upw-2*upwp) + wpf1.*(mean_u-mean_up) + upf1.*(mean_w-mean_wp) )/St;
ddt_vpwp   = ( mean_f1.*(wpv+vpw-2*vpwp) + wpf1.*(mean_v-mean_vp) + vpf1.*(mean_w-mean_wp) )/St;

% Terms upTp
ddt_upTp = ( mean_f1.*(Tpu-upTp) + Tpf1.*(mean_u-mean_up) )/St ...
    +(2*cr/(3*Pr*St))*( mean_f2.*(upT-upTp) + upf2.*(mean_T-mean_Tp) );
ddt_vpTp = ( mean_f1.*(Tpv-vpTp) + Tpf1.*(mean_v-mean_vp) )/St ...
    +(2*cr/(3*Pr*St))*( mean_f2.*(vpT-vpTp) + vpf2.*(mean_T-mean_Tp) );
ddt_wpTp = ( mean_f1.*(Tpw-wpTp) + Tpf1.*(mean_w-mean_wp) )/St ...
    +(2*cr/(3*Pr*St))*( mean_f2.*(wpT-wpTp) + wpf2.*(mean_T-mean_Tp) );

% Term cm2_Tp
ddt_cm2_Tp = (4*cr/(3*Pr*St))*( mean_f2.*(TpT-cm2_Tp) + Tpf2.*(mean_T-mean_Tp) );

% Terms xpalpha1
ddt_xpalpha1 = upalpha1;
ddt_ypalpha1 = vpalpha1;
ddt_zpalpha1 = wpalpha1;

% Terms upalpha1
ddt_upalpha1 = ( mean_f1.*(ualpha1-upalpha1) +f1alpha1.*(mean_u-mean_up) )/St;
ddt_vpalpha1 = ( mean_f1.*(valpha1-vpalpha1) +f1alpha1.*(mean_v-mean_vp) )/St;
ddt_wpalpha1 = ( mean_f1.*(walpha1-wpalpha1) +f1alpha1.*(mean_w-mean_wp) )/St;

% Terms xpalpha2
ddt_xpalpha2 = upalpha2;
ddt_ypalpha2 = vpalpha2;
ddt_zpalpha2 = wpalpha2;

% Terms upalpha12
ddt_upalpha2 = ( mean_f1.*(ualpha2-upalpha2) +f1alpha2.*(mean_u-mean_up) )/St;
ddt_vpalpha2 = ( mean_f1.*(valpha2-vpalpha2) +f1alpha2.*(mean_v-mean_vp) )/St;
ddt_wpalpha2 = ( mean_f1.*(walpha2-wpalpha2) +f1alpha2.*(mean_w-mean_wp) )/St;

% Term Tpalpha1
ddt_Tpalpha1 = (2*cr/(3*Pr*St))*( mean_f2.*(Talpha1-Tpalpha1) +f1alpha1.*(mean_T-mean_Tp) );

% Term Tpalpha2
ddt_Tpalpha2 = (2*cr/(3*Pr*St))*( mean_f2.*(Talpha2-Tpalpha2) +f1alpha2.*(mean_T-mean_Tp) );



% Outputs -----------------------------------------------------------------
out = [ ...
    ddt_mean_xp ; ddt_mean_yp ; ddt_mean_zp ; ...
    ddt_mean_up ; ddt_mean_vp ; ddt_mean_wp ; ...
    ddt_mean_Tp ; ...
    
    ddt_cm2_xp  ; ddt_cm2_yp  ; ddt_cm2_zp  ; ...
    ddt_xpyp    ; ddt_xpzp    ; ddt_ypzp    ; ...
    
    ddt_xpup    ; ddt_xpvp    ; ddt_xpwp    ; ...
    ddt_ypup    ; ddt_ypvp    ; ddt_ypwp    ; ...
    ddt_zpup    ; ddt_zpvp    ; ddt_zpwp    ; ...
    
    ddt_xpTp    ; ddt_ypTp    ; ddt_zpTp    ; ...
    
    ddt_cm2_up  ; ddt_cm2_vp  ; ddt_cm2_wp  ; ...
    ddt_upvp    ; ddt_upwp    ; ddt_vpwp    ; ...
    
    ddt_upTp    ; ddt_vpTp    ; ddt_wpTp    ; ...
    
    ddt_cm2_Tp  ; ...
    
    ddt_xpalpha1; ddt_ypalpha1; ddt_zpalpha1; ...
    ddt_upalpha1; ddt_vpalpha1; ddt_wpalpha1; ...
    ddt_xpalpha2; ddt_ypalpha2; ddt_zpalpha2; ...
    ddt_upalpha2; ddt_vpalpha2; ddt_wpalpha2; ...
    ddt_Tpalpha1; ddt_Tpalpha2 ...
    ];



end




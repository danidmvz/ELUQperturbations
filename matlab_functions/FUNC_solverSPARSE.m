%% SOLVER FOR 3D MoM ===================================================
% Note that here a is actually alpha, not the relative velocity a=u-u_p
% and the same for b, b=beta here, not the relative temperature b=T-T_p

function out = FUNC_solverSPARSE(inp)

% profile on
disp('FUNC_solverSPARSE')

% Obtaining inputs --------------------------------------------------------
% parameters
tLim       = inp.tLim;
nt         = inp.nt;
flowType   = inp.flowType;
TfieldType = inp.TfieldType;
f1Type     = inp.f1Type;
f2Type     = inp.f2Type;
timeMethod = inp.timeMethod;
Re_inf     = inp.Re_inf;
M_inf      = inp.M_inf;
dp         = inp.dp;
Pr         = inp.Pr;
cr         = inp.cr;
St         = inp.St;
nmp        = inp.nmp;
mean_a     = inp.mean_a;
mean_b     = inp.mean_b;
cm2_a      = inp.cm2_a;
ab         = inp.ab;

% Variables ---------------------------------------------------------------
% Means
mean_xp0 = inp.mean_xp0;
mean_yp0 = inp.mean_yp0;
mean_zp0 = inp.mean_zp0;
mean_up0 = inp.mean_up0;
mean_vp0 = inp.mean_vp0;
mean_wp0 = inp.mean_wp0;
mean_Tp0 = inp.mean_Tp0;

% xp with xp
cm2_xp0  = inp.cm2_xp0;
cm2_yp0  = inp.cm2_yp0;
cm2_zp0  = inp.cm2_zp0;
xpyp0    = inp.xpyp0;
xpzp0    = inp.xpzp0;
ypzp0    = inp.ypzp0;

% xp with up
xpup0    = inp.xpup0;
xpvp0    = inp.xpvp0;
xpwp0    = inp.xpwp0;
ypup0    = inp.ypup0;
ypvp0    = inp.ypvp0;
ypwp0    = inp.ypwp0;
zpup0    = inp.zpup0;
zpvp0    = inp.zpvp0;
zpwp0    = inp.zpwp0;

% xp with Tp
xpTp0    = inp.xpTp0;
ypTp0    = inp.ypTp0;
zpTp0    = inp.zpTp0;

% up with up
cm2_up0  = inp.cm2_up0;
cm2_vp0  = inp.cm2_vp0;
cm2_wp0  = inp.cm2_wp0;
upvp0    = inp.upvp0;
upwp0    = inp.upwp0;
vpwp0    = inp.vpwp0;

% up with Tp
upTp0    = inp.upTp0;
vpTp0    = inp.vpTp0;
wpTp0    = inp.wpTp0;

% Tp with Tp
cm2_Tp0  = inp.cm2_Tp0;

% a with xp
xpa0 = inp.xpa0;
ypa0 = inp.ypa0;
zpa0 = inp.zpa0;

% a with up
upa0 = inp.upa0;
vpa0 = inp.vpa0;
wpa0 = inp.wpa0;

% a with xp
xpb0 = inp.xpb0;
ypb0 = inp.ypb0;
zpb0 = inp.zpb0;

% a with up
upb0 = inp.upb0;
vpb0 = inp.vpb0;
wpb0 = inp.wpb0;

% a with Tp
Tpa0 = inp.Tpa0;

% b with Tp
Tpb0 = inp.Tpb0;



% Small calculations ------------------------------------------------------
nt_2save = nt(1);
nt_2iter = nt(2);
Nt       = (nt(1)-1)*(nt(2)-1)+1;
dt       = (tLim(2)-tLim(1))/(Nt-1);
t        = linspace(tLim(1),tLim(2),nt_2save)';

% Initializating ----------------------------------------------------------
time     = 0;
% Outputs -----------------------------------------------------------------
yyp = [ ...
    mean_xp0 ; mean_yp0 ; mean_zp0; ...
    mean_up0 ; mean_vp0 ; mean_wp0; ...
    mean_Tp0 ; ...
    
    cm2_xp0  ; cm2_yp0  ; cm2_zp0 ; ...
    xpyp0    ; xpzp0    ; ypzp0   ; ...
    
    xpup0    ; xpvp0    ; xpwp0   ; ...
    ypup0    ; ypvp0    ; ypwp0   ; ...
    zpup0    ; zpvp0    ; zpwp0   ; ...
    
    xpTp0    ; ypTp0    ; zpTp0   ; ...
    
    cm2_up0  ; cm2_vp0  ; cm2_wp0 ; ...
    upvp0    ; upwp0    ; vpwp0   ; ...
    
    upTp0    ; vpTp0    ; wpTp0   ; ...
    
    cm2_Tp0 ; ...
    xpa0  ; ypa0  ; zpa0 ; ...
    upa0  ; vpa0  ; wpa0 ;...
    xpb0  ; ypb0  ; zpb0 ; ...
    upb0  ; vpb0  ; wpb0 ;...
    Tpa0  ; Tpb0 ...
    ];

y   = yyp(:,1);


% Loop in time ------------------------------------------------------------
for j=1:nt_2save-1            % loop used to save in time
    for dummyVar=1:nt_2iter-1 % loop that iterates in every time step to save
        switch timeMethod
            case 1 % EE
                y  = y+dt*FUNC_rhsSPARSE(time,y,mean_a,mean_b,cm2_a,ab,flowType,TfieldType,f1Type,f2Type,Re_inf,M_inf,St,Pr,cr,dp,nmp);
            case 2 % RK2
                disp('not coded')
            case 3 % RK3 TVD
                fprintf('------------ RK3 TVD t=%.4f \n',time);
                k1 = FUNC_rhsSPARSE(time       ,                        y,mean_a,mean_b,cm2_a,ab,flowType,TfieldType,f1Type,f2Type,Re_inf,M_inf,St,Pr,cr,dp,nmp);
                k2 = FUNC_rhsSPARSE(time+dt    ,                  y+dt*k1,mean_a,mean_b,cm2_a,ab,flowType,TfieldType,f1Type,f2Type,Re_inf,M_inf,St,Pr,cr,dp,nmp);
                k3 = FUNC_rhsSPARSE(time+3/4*dt,y+(1/4)*dt*k1+(1/4)*dt*k2,mean_a,mean_b,cm2_a,ab,flowType,TfieldType,f1Type,f2Type,Re_inf,M_inf,St,Pr,cr,dp,nmp);
                y  = y+dt*((1/6)*k1+(1/6)*k2+(2/3)*k3);
            case 4 % RK4
                k1 = FUNC_rhsSPARSE(time       ,y           ,mean_a,mean_b,cm2_a,ab,flowType,TfieldType,f1Type,f2Type,Re_inf,M_inf,St,Pr,cr,dp,nmp);
                k2 = FUNC_rhsSPARSE(time+0.5*dt,y+0.5*dt.*k1,mean_a,mean_b,cm2_a,ab,flowType,TfieldType,f1Type,f2Type,Re_inf,M_inf,St,Pr,cr,dp,nmp);
                k3 = FUNC_rhsSPARSE(time+0.5*dt,y+0.5*dt.*k2,mean_a,mean_b,cm2_a,ab,flowType,TfieldType,f1Type,f2Type,Re_inf,M_inf,St,Pr,cr,dp,nmp);
                k4 = FUNC_rhsSPARSE(time+dt    ,y+    dt.*k3,mean_a,mean_b,cm2_a,ab,flowType,TfieldType,f1Type,f2Type,Re_inf,M_inf,St,Pr,cr,dp,nmp);
                y  = y+dt*(k1+2*k2+2*k3+k4)/6;
            otherwise
                disp('Error in timeMethod')
        end
        time = time+dt;
    end
    yyp(:,j+1) = y; % Saving after iterating without saving
    %     time
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
out.mean_xp = yyp( 0*nmp+1: 1*nmp,:)';
out.mean_yp = yyp( 1*nmp+1: 2*nmp,:)';
out.mean_zp = yyp( 2*nmp+1: 3*nmp,:)';
out.mean_up = yyp( 3*nmp+1: 4*nmp,:)';
out.mean_vp = yyp( 4*nmp+1: 5*nmp,:)';
out.mean_wp = yyp( 5*nmp+1: 6*nmp,:)';
out.mean_Tp = yyp( 6*nmp+1: 7*nmp,:)';

% xp with xp
out.cm2_xp  = yyp( 7*nmp+1: 8*nmp,:)';
out.cm2_yp  = yyp( 8*nmp+1: 9*nmp,:)';
out.cm2_zp  = yyp( 9*nmp+1:10*nmp,:)';
out.xpyp    = yyp(10*nmp+1:11*nmp,:)';
out.xpzp    = yyp(11*nmp+1:12*nmp,:)';
out.ypzp    = yyp(12*nmp+1:13*nmp,:)';

% xp with up
out.xpup    = yyp(13*nmp+1:14*nmp,:)';
out.xpvp    = yyp(14*nmp+1:15*nmp,:)';
out.xpwp    = yyp(15*nmp+1:16*nmp,:)';
out.ypup    = yyp(16*nmp+1:17*nmp,:)';
out.ypvp    = yyp(17*nmp+1:18*nmp,:)';
out.ypwp    = yyp(18*nmp+1:19*nmp,:)';
out.zpup    = yyp(19*nmp+1:20*nmp,:)';
out.zpvp    = yyp(20*nmp+1:21*nmp,:)';
out.zpwp    = yyp(21*nmp+1:22*nmp,:)';

% xp with Tp
out.xpTp    = yyp(22*nmp+1:23*nmp,:)';
out.ypTp    = yyp(23*nmp+1:24*nmp,:)';
out.zpTp    = yyp(24*nmp+1:25*nmp,:)';

% up with up
out.cm2_up  = yyp(25*nmp+1:26*nmp,:)';
out.cm2_vp  = yyp(26*nmp+1:27*nmp,:)';
out.cm2_wp  = yyp(27*nmp+1:28*nmp,:)';
out.upvp    = yyp(28*nmp+1:29*nmp,:)';
out.upwp    = yyp(29*nmp+1:30*nmp,:)';
out.vpwp    = yyp(30*nmp+1:31*nmp,:)';

% up with Tp
out.upTp    = yyp(31*nmp+1:32*nmp,:)';
out.vpTp    = yyp(32*nmp+1:33*nmp,:)';
out.wpTp    = yyp(33*nmp+1:34*nmp,:)';

% Tp with Tp
out.cm2_Tp  = yyp(34*nmp+1:35*nmp,:)';

% a with xp
out.xpa = yyp(35*nmp+1:36*nmp,:)';
out.ypa = yyp(36*nmp+1:37*nmp,:)';
out.zpa = yyp(37*nmp+1:38*nmp,:)';

% a with up
out.upa = yyp(38*nmp+1:39*nmp,:)';
out.vpa = yyp(39*nmp+1:40*nmp,:)';
out.wpa = yyp(40*nmp+1:41*nmp,:)';

% b with xp
out.xpb = yyp(41*nmp+1:42*nmp,:)';
out.ypb = yyp(42*nmp+1:43*nmp,:)';
out.zpb = yyp(43*nmp+1:44*nmp,:)';

% b with up
out.upb = yyp(44*nmp+1:45*nmp,:)';
out.vpb = yyp(45*nmp+1:46*nmp,:)';
out.wpb = yyp(46*nmp+1:47*nmp,:)';

% a with Tp
out.Tpa = yyp(47*nmp+1:48*nmp,:)';

% b with Tp
out.Tpb = yyp(48*nmp+1:49*nmp,:)';

% Time
out.t       = t;

% profile viewer

end
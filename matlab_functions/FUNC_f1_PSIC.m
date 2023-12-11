%% FUNCTION FOR DRAG CORRECTION f1 ========================================

function f1 = FUNC_f1_PSIC(u,v,w,up,vp,wp,f1Type,dataf1)

% Constructing Urel with units based on Uref ------------------------------
ax = u-up;
ay = v-vp;
az = w-wp;

switch f1Type
    
    case 1 % f1=1 (Stokes) ------------------------------------------------
        
        f1 = 1;
        
    case 2 % f1 = 1+0.15*Rep^0.687 ----------------------------------------
        
        modU  = sqrt(ax.^2 + ay.^2 + az.^2);
        Reinf = dataf1.Re_inf;
        dp    = dataf1.dp;
        Reinf = Re_inf;
        m     = 0.15;
        n     = 0.687;
        f1    = 1+m.*((ax.^2+ay.^2+az.^2).^(1/2).*dp.*Reinf).^n;
        
        
    case 3 % f1=1 (Stokes) random -----------------------------------------
        
        f1 = dataf1.a1(:);
        
    case 4 % f1 = 1+0.15*Rep^0.687 random ---------------------------------
        
        modU  = sqrt(ax.^2 + ay.^2 + az.^2);
        Reinf = dataf1.Re_inf;
        dp    = dataf1.dp;
        Reinf = Re_inf;
        m     = 0.15;
        n     = 0.687;
        f1    = dataf1.a1(:) .* ( 1+m.*((ax.^2+ay.^2+az.^2).^(1/2).*dp.*Reinf).^n ) ;
        
    otherwise
        disp('Error in f1Type');
end


end






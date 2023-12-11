%% FUNCTION FOR DRAG CORRECTION f1 ========================================

function out = FUNC_f1_SPARSE(u,v,w,up,vp,wp,f1Type,dataf1)

% Constructing Urel with units based on Uref ------------------------------
ax = u-up;
ay = v-vp;
az = w-wp;

switch f1Type
    
    case 1 % f1=1 (Stokes) ------------------------------------------------
        
        out.f1         = 1;
        out.df1dax     = 0;
        out.d2f1dax    = 0;
        out.d2f1daxday = 0;
        out.d2f1daxdaz = 0;
        
        out.df1day     = 0;
        out.d2f1day    = 0;
        out.d2f1daydaz = 0;
        
        out.df1daz     = 0;
        out.d2f1daz    = 0;
        
    case 2 % f1 = 1+0.15*Rep^0.687 ----------------------------------------
        
        modU  = sqrt(ax.^2 + ay.^2 + az.^2);
        Reinf = dataf1.Re_inf;
        dp    = dataf1.dp;
        Rep   = Reinf*modU*dp;
        Reinf = Re_inf;
        m     = 0.15;
        n     = 0.687;
        Recut = 0.1;
        out.f1         = 1+m.*((ax.^2+ay.^2+az.^2).^(1/2).*dp.*Reinf).^n;
        
        out.df1dax     = 0.*(Rep<Recut) + (Rep>=Recut)*(ax.*(ax.^2+ay.^2+az.^2).^(-1).*m.*n.*((ax.^2+ay.^2+az.^2).^(1/2).* ...
            dp.*Reinf).^n);
        out.d2f1dax    = 0.*(Rep<Recut) + (Rep>=Recut)*((ax.^2+ay.^2+az.^2).^(-2).*m.*(ay.^2+az.^2+ax.^2.*((-1)+n)).*n.*(( ...
            ax.^2+ay.^2+az.^2).^(1/2).*dp.*Reinf).^n);
        out.d2f1daxday = 0.*(Rep<Recut) + (Rep>=Recut)*(ax.*ay.*(ax.^2+ay.^2+az.^2).^(-2).*m.*((-2)+n).*n.*((ax.^2+ay.^2+ ...
            az.^2).^(1/2).*dp.*Reinf).^n);
        out.d2f1daxdaz = 0.*(Rep<Recut) + (Rep>=Recut)*(ax.*az.*(ax.^2+ay.^2+az.^2).^(-2).*m.*((-2)+n).*n.*((ax.^2+ay.^2+ ...
            az.^2).^(1/2).*dp.*Reinf).^n);
        
        out.df1day     = 0.*(Rep<Recut) + (Rep>=Recut)*(ay.*(ax.^2+ay.^2+az.^2).^(-1).*m.*n.*((ax.^2+ay.^2+az.^2).^(1/2).* ...
            dp.*Reinf).^n);
        out.d2f1day    = 0.*(Rep<Recut) + (Rep>=Recut)*((ax.^2+ay.^2+az.^2).^(-2).*m.*(ax.^2+az.^2+ay.^2.*((-1)+n)).*n.*(( ...
            ax.^2+ay.^2+az.^2).^(1/2).*dp.*Reinf).^n);
        out.d2f1daydaz = 0.*(Rep<Recut) + (Rep>=Recut)*(ay.*az.*(ax.^2+ay.^2+az.^2).^(-2).*m.*((-2)+n).*n.*((ax.^2+ay.^2+ ...
            az.^2).^(1/2).*dp.*Reinf).^n);
        
        out.df1daz  = 0.*(Rep<Recut) + (Rep>=Recut)*(az.*(ax.^2+ay.^2+az.^2).^(-1).*m.*n.*((ax.^2+ay.^2+az.^2).^(1/2).* ...
            dp.*Reinf).^n);
        out.d2f1daz = 0.*(Rep<Recut) + (Rep>=Recut)*((ax.^2+ay.^2+az.^2).^(-2).*m.*(ax.^2+ay.^2+az.^2.*((-1)+n)).*n.*(( ...
            ax.^2+ay.^2+az.^2).^(1/2).*dp.*Reinf).^n);
                   
    case 3 % f1=1 (Stokes) random -----------------------------------------
        
        out.f1         = dataf1.a1(:);
        out.df1dax     = 0;
        out.d2f1dax    = 0;
        out.d2f1daxday = 0;
        out.d2f1daxdaz = 0;
        
        out.df1day     = 0;
        out.d2f1day    = 0;
        out.d2f1daydaz = 0;
        
        out.df1daz     = 0;
        out.d2f1daz    = 0;
        
        
        
        
    otherwise
        disp('Error in f1Type');
end


end






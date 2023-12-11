%% FUNCTION FOR DRAG CORRECTION f1 ========================================

function out = FUNC_f2(time,u,v,w,up,vp,wp,f2Type,Re_inf,Pr,dp)

% Constructing Urel with units based on Uref ------------------------------
Urelx = u-up;
Urely = v-vp;
Urelz = w-wp;

% if Urelx == 0 & Urely == 0 & Urelz == 0
%     Urelx = eps;
%     Urely = eps;
%     Urelz = eps;
% end

switch f2Type
    case 1 % f2=1 (Laminar) -----------------------------------------------
        out.f2      = 1;
        out.df2dax  = 0;
        out.d2f2dax = 0;
        out.df2day  = 0;
        out.d2f2day = 0;
        out.df2daz  = 0;
        out.d2f2daz = 0;
        
    case 2 % Nu = 2+Re_p^0.5Pr^0.33; Nu=2*f2 ------------------------------
        
        Reinf       = Re_inf;
        Rho         = 1;
        
        out.f2      = 1+0.3E0.*Pr.^(1/3).*(dp.*Reinf.*Rho.*(Urelx.^2+Urely.^2+Urelz.^2) ...
            .^(1/2)).^(1/2);
        
        out.df2dax  = 0.15E0.*Pr.^(1/3).*Urelx.*(dp.*Reinf.*Rho.*(Urelx.^2+Urely.^2+ ...
            Urelz.^2).^(1/2)).^(1/2).*(Urelx.^2+Urely.^2+Urelz.^2).^(-1);
        
        out.d2f2dax = (-0.225E0).*Pr.^(1/3).*Urelx.^2.*(dp.*Reinf.*Rho.*(Urelx.^2+ ...
            Urely.^2+Urelz.^2).^(1/2)).^(1/2).*(Urelx.^2+Urely.^2+Urelz.^2).^( ...
            -2)+0.15E0.*Pr.^(1/3).*(dp.*Reinf.*Rho.*(Urelx.^2+Urely.^2+ ...
            Urelz.^2).^(1/2)).^(1/2).*(Urelx.^2+Urely.^2+Urelz.^2).^(-1);
        
        out.df2day  = 0.15E0.*Pr.^(1/3).*Urely.*(dp.*Reinf.*Rho.*(Urelx.^2+Urely.^2+ ...
            Urelz.^2).^(1/2)).^(1/2).*(Urelx.^2+Urely.^2+Urelz.^2).^(-1);
        
        out.d2f2day = (-0.225E0).*Pr.^(1/3).*Urely.^2.*(dp.*Reinf.*Rho.*(Urelx.^2+ ...
            Urely.^2+Urelz.^2).^(1/2)).^(1/2).*(Urelx.^2+Urely.^2+Urelz.^2).^( ...
            -2)+0.15E0.*Pr.^(1/3).*(dp.*Reinf.*Rho.*(Urelx.^2+Urely.^2+ ...
            Urelz.^2).^(1/2)).^(1/2).*(Urelx.^2+Urely.^2+Urelz.^2).^(-1);
        
        out.df2daz  = 0.15E0.*Pr.^(1/3).*Urelz.*(dp.*Reinf.*Rho.*(Urelx.^2+Urely.^2+ ...
            Urelz.^2).^(1/2)).^(1/2).*(Urelx.^2+Urely.^2+Urelz.^2).^(-1);
        
        out.d2f2daz = (-0.225E0).*Pr.^(1/3).*Urelz.^2.*(dp.*Reinf.*Rho.*(Urelx.^2+ ...
            Urely.^2+Urelz.^2).^(1/2)).^(1/2).*(Urelx.^2+Urely.^2+Urelz.^2).^( ...
            -2)+0.15E0.*Pr.^(1/3).*(dp.*Reinf.*Rho.*(Urelx.^2+Urely.^2+ ...
            Urelz.^2).^(1/2)).^(1/2).*(Urelx.^2+Urely.^2+Urelz.^2).^(-1);
        
        
%         if Urelx == 0 & Urely == 0 & Urelz == 0
%             out.df2dax  = 0;
%             out.d2f2dax = 0;
%             out.df2day  = 0;
%             out.d2f2day = 0;
%             out.df2daz  = 0;
%             out.d2f2daz = 0;
%         end
        
        if sqrt(Urelx.^2+Urely.^2+Urelz.^2)*Re_inf*dp<0.001
            out.f1      = 1;
            out.df2dax  = 0;
            out.d2f2dax = 0;
            out.df2day  = 0;
            out.d2f2day = 0;
            out.df2daz  = 0;
            out.d2f2daz = 0;
        end
        
%         Rep = Re_inf*sqrt(Urelx.^2+Urely.^2+Urelz.^2)*dp
        
        
    otherwise
        disp('Error in f1Type');
end


end






%% FLOW FUNCTION ==========================================================

function out = FUNC_flow(t,x,y,z,flowType,dataFlow)

np = length(x);

switch flowType
    
    case 0 % Uniform Flow -------------------------------------------------
        
        out.u(1:np,1) = ones(np,1);
        out.v(1:np,1) = zeros(np,1);
        out.w(1:np,1) = zeros(np,1);
        
        out.dudx(1:np,1)  = zeros(np,1);
        out.d2udx(1:np,1) = zeros(np,1);
        out.d2udxdy(1:np,1) = zeros(np,1);
        out.dudy(1:np,1)  = zeros(np,1);
        out.d2udy(1:np,1) = zeros(np,1);
        out.dudz(1:np,1)  = zeros(np,1);
        out.d2udz(1:np,1) = zeros(np,1);
        
        out.dvdx(1:np,1)  = zeros(np,1);
        out.d2vdx(1:np,1) = zeros(np,1);
        out.d2vdxdy(1:np,1) = zeros(np,1);
        out.dvdy(1:np,1)  = zeros(np,1);
        out.d2vdy(1:np,1) = zeros(np,1);
        out.dvdz(1:np,1)  = zeros(np,1);
        out.d2vdz(1:np,1) = zeros(np,1);
        
        out.dwdx(1:np,1)  = zeros(np,1);
        out.d2wdx(1:np,1) = zeros(np,1);
        out.dwdy(1:np,1)  = zeros(np,1);
        out.d2wdy(1:np,1) = zeros(np,1);
        out.dwdz(1:np,1)  = zeros(np,1);
        out.d2wdz(1:np,1) = zeros(np,1);
        
    case 3 % Double Gyre Flow ---------------------------------------------
        
        epsilon = dataFlow.e; % 0.1;
        w       = dataFlow.w;% 2*pi/10;
        A       = dataFlow.A;% 0.1;
        
        out.Psi = A.*sin(pi.*y).*sin(pi.*(epsilon.*x.^2.*sin(t.*w)+x.*(1+(-2).* ...
            epsilon.*sin(t.*w))));
        
        out.u(1:np,1)     = (-1).*A.*pi.*cos(pi.*y).*sin(pi.*x.*(1+epsilon.*((-2)+x).*sin(t.* ...
            w)));
        out.dudx(1:np,1)    = (-1).*A.*pi.^2.*cos(pi.*y).*cos(pi.*x.*(1+epsilon.*((-2)+x).*sin( ...
            t.*w))).*(1+2.*epsilon.*((-1)+x).*sin(t.*w));
        out.dudy(1:np,1)    = A.*pi.^2.*sin(pi.*y).*sin(pi.*x.*(1+epsilon.*((-2)+x).*sin(t.*w)));
        out.d2udx(1:np,1)   = A.*pi.^2.*cos(pi.*y).*((-2).*epsilon.*cos(pi.*x.*(1+epsilon.*((-2) ...
            +x).*sin(t.*w))).*sin(t.*w)+pi.*(1+2.*epsilon.*((-1)+x).*sin(t.*w) ...
            ).^2.*sin(pi.*x.*(1+epsilon.*((-2)+x).*sin(t.*w))));
        out.d2udy(1:np,1)   = A.*pi.^3.*cos(pi.*y).*sin(pi.*x.*(1+epsilon.*((-2)+x).*sin(t.*w)));
        out.d2udxdy(1:np,1) = A.*pi.^3.*cos(pi.*x.*(1+epsilon.*((-2)+x).*sin(t.*w))).*(1+2.* ...
            epsilon.*((-1)+x).*sin(t.*w)).*sin(pi.*y);
        
        out.v(1:np,1)       = A.*pi.*cos(pi.*x.*(1+epsilon.*((-2)+x).*sin(t.*w))).*(1+2.* ...
            epsilon.*((-1)+x).*sin(t.*w)).*sin(pi.*y);
        out.dvdx(1:np,1)    = A.*pi.*sin(pi.*y).*(2.*epsilon.*cos(pi.*x.*(1+epsilon.*((-2)+x).* ...
            sin(t.*w))).*sin(t.*w)+(-1).*pi.*(1+2.*epsilon.*((-1)+x).*sin(t.* ...
            w)).^2.*sin(pi.*x.*(1+epsilon.*((-2)+x).*sin(t.*w))));
        out.dvdy(1:np,1)    = A.*pi.^2.*cos(pi.*y).*cos(pi.*x.*(1+epsilon.*((-2)+x).*sin(t.*w))) ...
            .*(1+2.*epsilon.*((-1)+x).*sin(t.*w));
        out.d2vdx(1:np,1)   = A.*pi.^2.*(1+2.*epsilon.*((-1)+x).*sin(t.*w)).*sin(pi.*y).*((-1).* ...
            pi.*cos(pi.*x.*(1+epsilon.*((-2)+x).*sin(t.*w))).*(1+2.*epsilon.*( ...
            (-1)+x).*sin(t.*w)).^2+(-6).*epsilon.*sin(t.*w).*sin(pi.*x.*(1+ ...
            epsilon.*((-2)+x).*sin(t.*w))));
        out.d2vdy(1:np,1)   = (-1).*A.*pi.^3.*cos(pi.*x.*(1+epsilon.*((-2)+x).*sin(t.*w))).*(1+ ...
            2.*epsilon.*((-1)+x).*sin(t.*w)).*sin(pi.*y);
        out.d2vdxdy(1:np,1) = A.*pi.^2.*cos(pi.*y).*(2.*epsilon.*cos(pi.*x.*(1+epsilon.*((-2)+x) ...
            .*sin(t.*w))).*sin(t.*w)+(-1).*pi.*(1+2.*epsilon.*((-1)+x).*sin( ...
            t.*w)).^2.*sin(pi.*x.*(1+epsilon.*((-2)+x).*sin(t.*w))));
        
        
    case 4 % Sine flow --------------------------------------------------------
        
        k   = dataFlow.k;
        A   = dataFlow.A;
        w   = dataFlow.w;
        phi = dataFlow.phi;
        u0  = dataFlow.u0;
        out.u(1:np,1)     = u0+A*sin(phi-t*w+k*x);
        out.dudx(1:np,1)  = A.*k.*cos(phi+(-1).*t.*w+k.*x);
        out.d2udx(1:np,1) = (-1).*A.*k.^2.*sin(phi+(-1).*t.*w+k.*x);
        
        out.dudy(1:np,1)  = 0;
        out.d2udy(1:np,1) = 0;
        out.dudz(1:np,1)  = 0;
        out.d2udz(1:np,1) = 0;
        out.d2udxdy(1:np,1) = 0;
        out.d2udxdz(1:np,1) = 0;
        out.d2udydz(1:np,1) = 0;
        
        out.v(1:np,1)     = 0;
        out.dvdx(1:np,1)  = 0;
        out.d2vdx(1:np,1) = 0;
        out.dvdy(1:np,1)  = 0;
        out.d2vdy(1:np,1) = 0;
        out.dvdz(1:np,1)  = 0;
        out.d2vdz(1:np,1) = 0;
        out.d2vdxdy(1:np,1) = 0;
        out.d2vdxdz(1:np,1) = 0;
        out.d2vdydz(1:np,1) = 0;
        
        out.w(1:np,1)     = 0;
        out.dwdx(1:np,1)  = 0;
        out.d2wdx(1:np,1) = 0;
        out.dwdy(1:np,1)  = 0;
        out.d2wdy(1:np,1) = 0;
        out.dwdz(1:np,1)  = 0;
        out.d2wdz(1:np,1) = 0;
        out.d2wdxdy(1:np,1) = 0;
        out.d2wdxdz(1:np,1) = 0;
        out.d2wdydz(1:np,1) = 0;
end

end


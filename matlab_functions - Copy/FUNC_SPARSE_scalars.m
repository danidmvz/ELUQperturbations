%%

function S1 = FUNC_SPARSE_scalars(S1,alpha)


for j=1:S1.inp.nt(1)
    
    S1.cx(:,j)      = S1.mean_xp(:,j);
    S1.cy(:,j)      = S1.mean_yp(:,j);
    [S1.L(:,j),~,~] = curvature([S1.cx(:,j) S1.cy(:,j)]);
    
    alpha = S1.L(:,j);
    
    %     S1.dds_cov_ang(:,j)      = FUNC_derivative(S1.cov_ang(:,j),alpha);
    %     S1.dds_cov_ang(:,j)      = FUNC_derivative(S1.cov_ang(:,j),alpha);
    
    
    for i=1:length(alpha)
        S1.cov(:,:,i,j)     = [S1.cm2_xp(i,j) S1.xpyp(i,j) ; S1.xpyp(i,j) S1.cm2_yp(i,j)];
        S1.max_cov_eig(i,j) = max(abs(eig(S1.cov(:,:,i,j))));
        [S1.cov_V(:,:,j),~] = eig(S1.cov(:,:,i,j));
        S1.cov_ang(i,j)     = atan(S1.cov_V(2,1,j)/S1.cov_V(1,1,j));
        
        S1.ddt_cov(:,:,i,j)     = [S1.RHS.cm2_xp(i,j) S1.RHS.xpyp(i,j) ; S1.RHS.xpyp(i,j) S1.RHS.cm2_yp(i,j)];
        S1.ddt_max_cov_eig(i,j) = max(abs(eig(S1.ddt_cov(:,:,i,j))));
        
        S1.sigma(i,j)    = log(sqrt(S1.max_cov_eig(i,j)./S1.max_cov_eig(i,1)));
        S1.ddt_sigma(i,j)= log(sqrt(S1.ddt_max_cov_eig(i,j)./S1.ddt_max_cov_eig(i,1)));
        
        % Product of eigenvalues
        d                    = eig(S1.cov(:,:,i,j));
        S1.prod_cov_eig(i,j) = d(1)*d(2);
    end
    
    
    S1.dds_cov_ang(:,j)      = FUNC_derivative(S1.cov_ang(:,j),alpha);
    S1.dds_max_cov_eig(:,j)  = abs(FUNC_derivative(S1.max_cov_eig(:,j),alpha));
    S1.d2ds_max_cov_eig(:,j) = abs((FUNC_derivative(S1.dds_max_cov_eig(:,j),alpha)));
    
    S1.dds_ddt_max_cov_eig(:,j) = FUNC_derivative(S1.ddt_max_cov_eig(:,j),alpha);
    S1.dds_sigma(:,j)           = FUNC_derivative(S1.sigma(:,j),alpha);
    
    S1.dds_xpyp(:,j)  = abs(FUNC_derivative(S1.xpyp(:,j),alpha));
    
end
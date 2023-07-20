%% --- 
% Adjoint Transformation:
% - depends on current configuration of the manipulator
% --- 
A_ij = {};
for i=1:N_jnts
    for j=1:N_jnts
        % disp([i,j])
        if i>j
            G_ji = Lie.prod_SE3_from_SE3xk({exp_xi_theta_in_SE3_{j+1:i}});
            A_ij{i,j} = Lie.inv_Ad_SE3_from_SE3(G_ji);
        elseif i==j
            A_ij{i,j} = eye(6);
        else
            A_ij{i,j} = 0;
        end
    end
end

%% --- 
% Compute Moments and Inertia:
% --- 
% dtheta = symmatrix('dt', [N_jnts,1])
dtheta = zeros(N_jnts,1);

validate.set_validatorLevel("false")
M_ij_ = [];
C_ij_ = [];
for i=1:N_jnts
    for j=1:N_jnts
        M_ij = 0;
        for l=max(i,j):N_jnts
            % [Murray 4.29] Inertia Matrix:
            M_ij = M_ij ...
                + xi_R6_0_{i}' * A_ij{l,i}' * M_R6x6_spatial_{l} * A_ij{l,j} * xi_R6_0_{j};
        end
        M_ij_(i,j) = M_ij;
        % M_ij_{i,j} = M_ij;
        
        C_ij = 0;
        for k=1:N_jnts
            % [Murray 4.29] Coriolis Matrix:
            dM_dt = 0;
            if k == 1
                A1 = 0;
                A2 = 0;
            else
                A1 = A_ij{k-1,i};
                A2 = A_ij{k-1,j};
            end
            
            % [Murray 4.30] Jacobian of Moments:
            A1_xi_k = Lie.lie_bracket_from_R6xR6( A1 * xi_R6_0_{i}, xi_R6_0_{k} )';
            A2_xi_k = Lie.lie_bracket_from_R6xR6( A2 * xi_R6_0_{j}, xi_R6_0_{k} )';
            
            A1_xi_j = Lie.lie_bracket_from_R6xR6( A1 * xi_R6_0_{i}, xi_R6_0_{j} )';
            A2_xi_j = Lie.lie_bracket_from_R6xR6( A2 * xi_R6_0_{k}, xi_R6_0_{j} )';
            
            A1_xi_i = Lie.lie_bracket_from_R6xR6( A1 * xi_R6_0_{k}, xi_R6_0_{i} )';
            A2_xi_i = Lie.lie_bracket_from_R6xR6( A2 * xi_R6_0_{j}, xi_R6_0_{i} )';

            dMij_dtk = 0;
            dMik_dtj = 0;
            dMkj_dti = 0;
            for l=max(i,j):N_jnts
                dMij_dtk = dMij_dtk ...
                    + A1_xi_k' * A_ij{l,k}' * M_R6x6_spatial_{l} * A_ij{l,j} * xi_R6_0_{j} ...
                    + xi_R6_0_{i}' * A_ij{l,i}' * M_R6x6_spatial_{l} * A_ij{l,k} * A2_xi_k;
                
                dMik_dtj = dMik_dtj ...
                    + A1_xi_j' * A_ij{l,j}' * M_R6x6_spatial_{l} * A_ij{l,k} * xi_R6_0_{k} ...
                    + xi_R6_0_{i}' * A_ij{l,i}' * M_R6x6_spatial_{l} * A_ij{l,j} * A2_xi_j;

                dMkj_dti = dMkj_dti ...
                    + A1_xi_i' * A_ij{l,i}' * M_R6x6_spatial_{l} * A_ij{l,j} * xi_R6_0_{j} ...
                    + xi_R6_0_{k}' * A_ij{l,k}' * M_R6x6_spatial_{l} * A_ij{l,i} * A2_xi_i;
            end
            C_ij = C_ij + ( (dMij_dtk + dMik_dtj + dMkj_dti) * dtheta(k) );
        end
        % C_ij_{i,j} = 0.5 * C_ij; % cache as struct for symbolic
        C_ij_(i,j) = 0.5 * C_ij; % cache as matrix for numerical
    end
end

validate.set_validatorLevel("all")
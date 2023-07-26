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


%% % otheres
function [G_SE3_wam_spatial_0_, J_spatial_] = compute_Spatial(G_SE3_0_spatial, xi_R6_0_, theta_)
    J_spatial_ = [];
    G_SE3_wam_spatial_0_ = {};
    
    % init:
    N_jnts = length(theta_);
    G_SE3_wam_spatial_0_{1} = G_SE3_0_spatial; % wam base location
    
    % Compute Spatial (Forward):
    for i=1:N_jnts
        % -> compute transformation from twist coordinates and angles applied:
        % |---- (more general):
        G_i_t_{i} = RBT.screw_SE3_from_xi_theta(xi_R6_0_{i}, theta_(i));
        % |---- (more efficient):
        % G_i_t_{i} = RBT.screw_SE3_from_axis_point_theta(w_R3_0_{i}, q_R3_0_{i}, theta_(i));
        
        % -> cumulate transformation:
        G_SE3_wam_spatial_0_{i+1} = G_SE3_wam_spatial_0_{i} * G_i_t_{i};
        
        % -> compute jacobian:
        J_spatial_i = Lie.Ad_SE3_from_SE3(G_SE3_wam_spatial_0_{i}) * xi_R6_0_{i};
        % output:
        J_spatial_ = [J_spatial_, J_spatial_i];
    end
end
function [G_SE3_wam_body_0_, J_body_] = compute_Body(G_SE3_tool_frame, xi_R6_0_, theta_)
    J_body_ = [];
    G_SE3_wam_body_0_ = {};
    
    % init:
    N_jnts = length(theta_);
    G_SE3_wam_body_0_{1} = G_SE3_tool_frame;
    
    %% Compute Body (Backward):
    for i=1:N_jnts
        % - compute backwards from tool to base:
        G_b_t_{i} = RBT.screw_SE3_from_xi_theta(xi_R6_0_{N_jnts-i+1}, theta_(N_jnts-i+1));
        % -> cumulate transformation:
        G_SE3_wam_body_0_{i+1} = G_SE3_wam_body_0_{i} * G_b_t_{i};
        
        % -> compute jacobian:
        J_body_i = Lie.inv_Ad_SE3_from_SE3(G_SE3_wam_body_0_{i}) * xi_R6_0_{N_jnts-i+1};
        % output:
        J_body_ = [J_body_i, J_body_];
    end
end
function G_SE_b_st_1_ = compute_Body_FK(G_SE3_0_st, exp_xi_theta_in_SE3_1_)
    % @G_SE3_0_st : zero/home configuration of the EE
    N_jnts = length(exp_xi_theta_in_SE3_1_);
    % Compute Body (Forward):
    % > g_st(t) = exp_xi_1(t_i) * exp_xi_2(t_i+1) * ... * exp_xi_n(t_n) * g_st(0)
    T = G_SE3_0_st;
    for i=1:N_jnts
        T = exp_xi_theta_in_SE3_1_{N_jnts-i+1} * T; % cummulate
        G_SE3_s_st_1_{i} = T;                       % cache
    end
    G_SE3_s_st_1_ = flip(G_SE3_s_st_1_); % flip the order
end
function [G_SE3_s_st_, J_6x6_s_st_] = compute_Spatial_from_SE3(G_SE3_0_spatial, xi_R6_0_, exp_xi_theta_in_SE3_)
    % init:
    N_jnts = length(exp_xi_theta_in_SE3_);
    J_6x6_s_st_ = [];                % joints starting from 1 ... N
    G_SE3_s_st_ = {G_SE3_0_spatial}; % frames starting from 0 ... N
    
    % Compute Spatial (Forward):
    for i=1:N_jnts
        % -> compute current jacobian:
        Ad_0_i = Lie.Ad_SE3_from_SE3(G_SE3_s_st_{i}); % frame i-1:0
        J_6x6_s_st_ =  Ad_0_i * xi_R6_0_{i};          % joint   i:1

        % -> cumulate current transformation:
        G_SE3_s_st_{i+1} = G_SE3_s_st_{i} * exp_xi_theta_in_SE3_{i}; % frame i:1
    end
end


%% %%
G_SE3_wam_spatial_0 = OpenChain.compute_Spatial_FK(exp_xi_theta_in_SE3_, G_SE3_0_{end})
% G_SE3_wam_spatial_0 = OpenChain.compute_Spatial_FK_from_theta(xi_R6_1_, JOINT_ANGLEG_SE3_0_{end})

% [G_SE3_wam_spatial_0_, J_spatial_] = OpenChain.compute_Spatial_from_SE3(SUMMIT_POSE_SE3, xi_R6_1_, exp_xi_theta_in_SE3_);
% [G_SE3_wam_body_0_, J_body_] = OpenChain.compute_Body_from_SE3(G_SE3_0_{end}, xi_R6_1_, exp_xi_theta_in_SE3_);

% -> alternative: compute directly
% [G_SE3_wam_spatial_, J_spatial_] = OpenChain.compute_Spatial(G_SE3_summit_wam, xi_R6_0_, JOINT_ANGLE);
% [G_SE3_wam_body_, J_body_] = OpenChain.compute_Body(G_SE3_0_{N_jnts}, xi_R6_0_, JOINT_ANGLE);

%% --- 

%
theta = zeros(N_jnts);
d_theta = zeros(N_jnts);
dd_theta = zeros(N_jnts);

g = [0; 0; 0];
Ftip = [0; 0; 0; 0; 0; 0]; % > Spatial force applied by the end-effector expressed in frame {n+1}
Mlist = cat(3, G_SE3_{:});
Glist = cat(3, M_R6x6_{:});

Slist = cat(2, xi_R6_0_{:});


Js = OpenChainMR.spatial_jacobian(Slist, theta);

% Tau = OpenChainMR.inverse_dynamics(theta, d_theta, dd_theta, g, Ftip, Mlist, Glist, Slist);
% c = OpenChainMR.vel_qualdratic_force(theta, d_theta, Mlist, Glist, Slist);
% M = OpenChainMR.mass_matrix(theta, Mlist, Glist, Slist)

% Jb = OpenChainMR.body_jacobian(Blist, theta)

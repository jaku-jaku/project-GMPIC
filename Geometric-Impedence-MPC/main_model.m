close all;
clear all;
clc;

%% [INIT] ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 
% --- 
addpath(genpath("utils"));
% --- 
helper.logtitle("Initialize")
% --- 
% - USER PARAMs:
SAVE_CONSOLE = false;
CLEAR_OUTPUT = true;
CLOSE_WINDOW = true; % pre-closing
AUTO_CLOSE   = false;

helper.createFolder("output/test", false);
helper.setLogLevel("all")
validate.set_validatorLevel("none")
% --- 
helper.endSection(AUTO_CLOSE);
%% Initial g(0)) ===== ===== ===== ===== ===== ===== =====:
DIR = helper.declareSection("test", "init_wam", SAVE_CONSOLE, CLEAR_OUTPUT, CLOSE_WINDOW);
% --- 

% [ i-->i+1 Links ]
N_jnts = length(common.WAM)-1;
for i = 1:N_jnts
    helper.loginfo(sprintf("> Indexing Links @ %d-->%d",i-1,i));
    % - obtain relative params:
    q_R3_i = common.WAM(i).tail_frame_transformation.q';
    w_R3_i = common.WAM(i).tail_frame_transformation.w';
    t_R_i  = common.WAM(i).tail_frame_transformation.t;
    w_joint_axis_i = common.WAM(i).tail_frame_revolute_joint.axis';

    % -> (axis,angle) to SO3:
    R_SO3_i = Lie.rodrigues_SO3_from_unit_R3xR(w_R3_i, t_R_i);

    JOINT_CAD_AXIS = [0;0;1]; % Along link-1 z-axis
    % - on SE3 operation:
    % - grab local frame:
    q_R4_i = [q_R3_i; 1]; % relative displacement
    G_SE3_i = [R_SO3_i, q_R3_i; zeros(1,3), 1]; % relative RBT
    % - convert to base frame:
    if i > 1
        q_R4_0_i = G_SE3_0_{i-1} * q_R4_i;
        G_SE3_0_i = G_SE3_0_{i-1} * G_SE3_i;
    else
        q_R4_0_i = eye(4) * q_R4_i;
        G_SE3_0_i = eye(4) * G_SE3_i;
    end
    % - obtain global form of local z-axis
    w_R3_0_i = G_SE3_0_i(1:3,1:3) * w_joint_axis_i; 
    
    % - twist:
    xi_R6_0_i = RBT.twist_R6_from_axis_point(w_R3_0_i,q_R4_0_i(1:3));
    helper.logdebug(helper.a2str("xi_R6_0_i",xi_R6_0_i));
    
    % - compute generalized inertia matrix:
    m_R_i = common.WAM(i).mass;
    mc_R3_i = common.WAM(i).tail_frame_mass_center';
    I_mc_R3x3_i = common.WAM(i).Tail_frame_MoI_at_mass; % inertia at mass center , aligned with output frame axis
    mc_R3_i_hat = Lie.hat_so3_from_R3(mc_R3_i);
    % - inertia at mass center , aligned with output frame axis:
    % M_mass_center_i = [
    %     m_R_i * eye(3)  , zeros(3,3);
    %     zeros(3,3)      , I_mc_R3x3_i;
    % ];
    % - translating inertia matrix from center-of-mass to the output frame:
    M_R6x6_{i} = [
        m_R_i * eye(3)       , -m_R_i*mc_R3_i_hat; ...
        m_R_i * mc_R3_i_hat  , I_mc_R3x3_i - m_R_i * mc_R3_i_hat^2
    ]; % [pg 288, ]

    % - adjoint (0)
    % Ad_inv_g0_s_l_{i} = Lie.Ad_SE3_from_SE3(RBT.inverse_SE3(G_SE3_0_i));
    Ad_inv_g0_s_l_{i} = Lie.inv_Ad_SE3_from_SE3(G_SE3_0_i); % equivalent

    % - Inertia of the ith link reflected into the base spatial frame:
    % [4.28, Murray]
    M_R6x6_spatial_{i} = Ad_inv_g0_s_l_{i}' * M_R6x6_{i} * Ad_inv_g0_s_l_{i}; 
    
    % (cache):
    w_R3_0_{i} = w_R3_0_i;
    q_R3_0_{i} = q_R4_0_i(1:3);
    xi_R6_0_{i} = xi_R6_0_i;
    G_SE3_0_{i} = G_SE3_0_i; % transformation
    helper.logdebug(helper.a2str("G_SE3_0_i",G_SE3_0_i));
end

% [ Define ]

%% Compute Kinematics) ===== ===== ===== ===== ===== ===== =====:
DIR = helper.declareSection("test", "compute", SAVE_CONSOLE, CLEAR_OUTPUT, CLOSE_WINDOW);
% --- 
% -> (axis,angle) to SO3:
R_SO3_summit = Lie.rodrigues_SO3_from_unit_R3xR([0;0;1], 0);
% (homogeneous):
G_SE3_summit_wam = [R_SO3_summit, common.dP_SUMMIT_WAM.'; zeros(1,3), 1]; % relative RBT

% Given Angles:
JOINT_ANGLE = ones(1,N_jnts) * pi/4;
% compute FK:
exp_xi_theta_in_SE3_ = OpenChain.batch_screw_SE3_from_twist_angle_R6xR(xi_R6_0_, JOINT_ANGLE);

% --- 
% spacial jacobian:
% body jacobian:
% --- 
[G_SE3_wam_spatial_0_, J_spatial_] = OpenChain.compute_Spatial_from_SE3(G_SE3_summit_wam, xi_R6_0_, exp_xi_theta_in_SE3_);
[G_SE3_wam_body_0_, J_body_] = OpenChain.compute_Body_from_SE3(G_SE3_0_{N_jnts}, xi_R6_0_, exp_xi_theta_in_SE3_);

% [G_SE3_wam_spatial_, J_spatial_] = OpenChain.compute_Spatial(G_SE3_summit_wam, xi_R6_0_, JOINT_ANGLE);
% [G_SE3_wam_body_, J_body_] = OpenChain.compute_Body(G_SE3_0_{N_jnts}, xi_R6_0_, JOINT_ANGLE);

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
dtheta = ones(N_jnts,1);
dtheta(end) = 0;

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

% A = blkdiag(xi_R6_0_{:});
% G = blkdiag(M_R6x6_spatial_{:});
% 



%% PLOT) ===== ===== ===== ===== ===== ===== =====:
helper.endSection(AUTO_CLOSE);
DIR = helper.declareSection("test", "plot_wam", SAVE_CONSOLE, CLEAR_OUTPUT, CLOSE_WINDOW);
helper.newFigure(-1);
% tiledlayout(2,1);
ax_1_1 = nexttile();

% plot summit:
utils.plot_Summit(common.SUMMIT_INIT_POSE,common.SUMMIT_INIT_POSE,'S',ax_1_1);

% plot wam base link:
utils.plot_link( ...
        G_SE3_summit_wam, ...
        eye(4), ...
        common.WAM(1).tail_frame_transformation.q', ...
        common.WAM(1).base_frame, ...
        'k', ax_1_1);

% plot rest links:
for i=1:N_jnts
    % -> plot:
    utils.plot_link( ...
        G_SE3_wam_spatial_0_{i+1}, ...
        G_SE3_0_{i}, ...
        common.WAM(i+1).tail_frame_transformation.q', ...
        common.WAM(i+1).base_frame, ...
        'k', ax_1_1);
end

grid on;
axis equal;
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')
% axis([-1 1 -1 1 0 1])
helper.saveFigure([400,600], DIR, "FK")

% --- 
helper.endSection(AUTO_CLOSE);
%% [EOF] ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 

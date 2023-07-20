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
validate.set_validatorLevel("all")
% --- 
helper.endSection(AUTO_CLOSE);
%% Initial g(0)) ===== ===== ===== ===== ===== ===== =====:
DIR = helper.declareSection("test", "init_wam", SAVE_CONSOLE, CLEAR_OUTPUT, CLOSE_WINDOW);
% --- 


% Summit ---> WAM base frame:
R_SO3_summit = Lie.rodrigues_SO3_from_unit_R3xR([0;0;1], 0); % -> (axis,angle) to SO3:
G_SE3_summit_wam = [R_SO3_summit, common.dP_SUMMIT_WAM.'; zeros(1,3), 1]; % relative RBT

% define base spatial frame for WAM:
WAM_Spatial_0 = G_SE3_summit_wam;
% WAM_Spatial_0 = eye(4);

% [ i-1-->i: i's Links with i's output frame]
N_links = length(common.WAM);
N_jnts = N_links; % we cosider (i-1)===>(i)'s joints/frames
valid_jnts_filter = eye(N_jnts); % <--- assume all joints active by default

G_SE3_0_{1} = WAM_Spatial_0; % <--- initial wam base frame

for i = 1:N_links
    helper.loginfo(sprintf("> Indexing Links @ %d-->%d",i-1,i));
    % - obtain relative params:
    q_R3_i = common.WAM(i).tail_frame_transformation.q';
    w_R3_i = common.WAM(i).tail_frame_transformation.w';
    t_R_i  = common.WAM(i).tail_frame_transformation.t;
    w_joint_axis_i = common.WAM(i).tail_frame_revolute_joint.axis';

    joint_name_i = common.WAM(i).tail_frame_revolute_joint.name;
    joint_limits_{i} = common.WAM(i).tail_frame_revolute_joint.limits;

    % - store joint info:
    if joint_name_i == "NotActive"
        valid_jnts_filter(i,i) = 0; % invalid angle inputs
    end

    % -> (axis,angle) to SO3:
    R_SO3_i = Lie.rodrigues_SO3_from_unit_R3xR(w_R3_i, t_R_i);

    JOINT_CAD_AXIS = [0;0;1]; % Along link-1 z-axis
    % - on SE3 operation:
    % - grab link tail frame transformation in local frame:
    q_R4_{i} = [q_R3_i; 1]; % relative displacement
    G_SE3_{i} = [R_SO3_i, q_R3_i; zeros(1,3), 1]; % relative RBT
    % - convert to link base frame:
    q_R4_0_i = G_SE3_0_{i} * q_R4_{i};
    G_SE3_0_i = G_SE3_0_{i} * G_SE3_{i};

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
    G_SE3_0_{i+1} = G_SE3_0_i; % transformation
    helper.logdebug(helper.a2str("G_SE3_0_i",G_SE3_0_i));
end

% [ Define ]

%% Compute Kinematics) ===== ===== ===== ===== ===== ===== =====:
DIR = helper.declareSection("test", "compute", SAVE_CONSOLE, CLEAR_OUTPUT, CLOSE_WINDOW);
% --- 

% Given Angles:
% JOINT_ANGLE = ones(1,N_jnts) * pi/4 * valid_jnts_filter;
JOINT_ANGLE =zeros(1,N_jnts);
% compute FK:
exp_xi_theta_in_SE3_ = OpenChain.batch_screw_SE3_from_twist_angle_R6xR(xi_R6_0_, JOINT_ANGLE);

% --- 
% spacial jacobian:
% body jacobian:
% --- 
SUMMIT_POSE_SE3 = eye(4) * common.SUMMIT_INIT_POSE; % summit stationary --> WAM stationary
[G_SE3_wam_spatial_0_, J_spatial_] = OpenChain.compute_Spatial_from_SE3(SUMMIT_POSE_SE3, xi_R6_0_, exp_xi_theta_in_SE3_);
[G_SE3_wam_body_0_, J_body_] = OpenChain.compute_Body_from_SE3(G_SE3_0_{end}, xi_R6_0_, exp_xi_theta_in_SE3_);

% [G_SE3_wam_spatial_, J_spatial_] = OpenChain.compute_Spatial(G_SE3_summit_wam, xi_R6_0_, JOINT_ANGLE);
% [G_SE3_wam_body_, J_body_] = OpenChain.compute_Body(G_SE3_0_{N_jnts}, xi_R6_0_, JOINT_ANGLE);

%% --- 

%
% A = blkdiag(xi_R6_0_{:});
% G = blkdiag(M_R6x6_spatial_{:});
% 
% Takes thetalist: n-vector of joint variables,
%       dthetalist: n-vector of joint rates,
%       ddthetalist: n-vector of joint accelerations,
%       g: Gravity vector g,
%       Ftip: Spatial force applied by the end-effector expressed in frame 
%             {n+1},
%       Mlist: List of link frames {i} relative to {i-1} at the home 
%              position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format
%              of a matrix with the screw axes as the columns.
% Returns taulist: The n-vector of required joint forces/torques.
% This function uses forward-backward Newton-Euler iterations to solve the 
% equation:
% taulist = Mlist(thetalist) * ddthetalist + c(thetalist, dthetalist) ...
%           + g(thetalist) + Jtr(thetalist) * Ftip
% Example Input (3 Link Robot):
theta = zeros(N_jnts);
d_theta = dtheta;
dd_theta = zeros(N_jnts);

g = [0; 0; 0];
Ftip = [0; 0; 0; 0; 0; 0]; % > Spatial force applied by the end-effector expressed in frame {n+1}
Mlist = cat(3, G_SE3_{:});
Glist = cat(3, M_R6x6_{:});

Slist = cat(2, xi_R6_0_{:});


Js = OpenChain.space_jacobian(Slist, theta)

% Tau = OpenChain.inverse_dynamics(theta, d_theta, dd_theta, g, Ftip, Mlist, Glist, Slist);
% c = OpenChain.vel_qualdratic_force(theta, d_theta, Mlist, Glist, Slist);
% M = OpenChain.mass_matrix(theta, Mlist, Glist, Slist)

% Jb = OpenChain.body_jacobian(Blist, theta)




%% PLOT) ===== ===== ===== ===== ===== ===== =====:
helper.endSection(AUTO_CLOSE);
DIR = helper.declareSection("test", "plot_wam", SAVE_CONSOLE, CLEAR_OUTPUT, CLOSE_WINDOW);
helper.newFigure(-1);
% tiledlayout(2,1);
ax_1_1 = nexttile();

% plot summit:
utils.plot_Summit(SUMMIT_POSE_SE3,SUMMIT_POSE_SE3,'S',ax_1_1);

% plot links:
for i=1:N_links
    % -> plot:
    if i == N_links 
        % plot the end-effector frame as well:
        utils.plot_link( ...
            G_SE3_0_{i}, ...    
            G_SE3_0_{i+1}, ...
            G_SE3_wam_spatial_0_{i}, ...
            G_SE3_wam_spatial_0_{i+1}, ...
            common.WAM(i).tail_frame_transformation.q', ...
            common.WAM(i).base_frame, ...
            common.WAM(i).tail_frame, ...
            'k', ax_1_1);
    else
        utils.plot_link( ...
            G_SE3_0_{i}, ...    
            [], ...
            G_SE3_wam_spatial_0_{i}, ...
            [], ...
            common.WAM(i).tail_frame_transformation.q', ...
            common.WAM(i).base_frame, ...
            [], ...
            'k', ax_1_1);
    end
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

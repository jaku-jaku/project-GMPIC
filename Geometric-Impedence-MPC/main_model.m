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
%% Tests) ===== ===== ===== ===== ===== ===== =====:
DIR = helper.declareSection("test", "init_wam", SAVE_CONSOLE, CLEAR_OUTPUT, CLOSE_WINDOW);
% --- 

% [ i-->i+1 Links ]
N_jnts = length(common.WAM)-1;
for i = 1:N_jnts
    helper.loginfo(sprintf("> Indexing Links @ %d-->%d",i-1,i));
    % - obtain relative params:
    q_R3_i = common.WAM(i).tip_frame_transformation.q';
    w_R3_i = common.WAM(i).tip_frame_transformation.w';
    t_R_i  = common.WAM(i).tip_frame_transformation.t;
    w_joint_axis_i = common.WAM(i).tip_frame_revolute_joint.axis';

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
    
    % (cache):
    w_R3_0_{i} = w_R3_0_i;
    q_R3_0_{i} = q_R4_0_i(1:3);
    xi_R6_0_{i} = xi_R6_0_i;
    G_SE3_0_{i} = G_SE3_0_i; % transformation
    helper.logdebug(helper.a2str("G_SE3_0_i",G_SE3_0_i));
end

% [ Define ]

% --- 
helper.endSection(AUTO_CLOSE);
%% Compute Kinematics) ===== ===== ===== ===== ===== ===== =====:
DIR = helper.declareSection("test", "plot_wam", SAVE_CONSOLE, CLEAR_OUTPUT, CLOSE_WINDOW);
% --- 
% -> (axis,angle) to SO3:
R_SO3_summit = Lie.rodrigues_SO3_from_unit_R3xR([0;0;1], 0);
% (homogeneous):
G_SE3_summit_wam = [R_SO3_summit, common.dP_SUMMIT_WAM.'; zeros(1,3), 1]; % relative RBT

% special position, without joint limits applied:
JOINT_ANGLE = zeros(N_jnts);

% Compute Spatial (Forward):
J_spatial = [];
G_SE3_wam_spatial_ = {};
G_prev = G_SE3_summit_wam;
for i=1:N_jnts
    theta = JOINT_ANGLE(i);
    % -> compute transformation from twist coordinates and angles applied:
    % |---- (more general):
    G_i_t_{i} = RBT.screw_SE3_from_xi_theta(xi_R6_0_{i}, theta);
    % |---- (more efficient):
    % G_i_t_{i} = RBT.screw_SE3_from_axis_point_theta(w_R3_0_{i}, q_R3_0_{i}, theta);
    
    % -> cumulate transformation:
    G_SE3_wam_spatial_i = G_prev * G_i_t_{i};
    G_prev = G_SE3_wam_spatial_i;
    % output:
    G_SE3_wam_spatial_{i} = G_SE3_wam_spatial_i;

    % -> compute jacobian:
    J_spatial_i = Lie.Ad_SE3_from_SE3(G_prev) * xi_R6_0_{i};
    % output:
    J_spatial = [J_spatial, J_spatial_i];
end

%% Compute Body (Backward):
J_body = [];
G_SE3_wam_body_ = {};
G_prev = G_SE3_0_{N_jnts};
for i=1:N_jnts
    G_b_t_{i} = RBT.screw_SE3_from_xi_theta(xi_R6_0_{N_jnts-i+1}, theta);
    % -> cumulate transformation:
    G_SE3_wam_body_i = G_prev * G_b_t_{i};
    G_prev = G_SE3_wam_body_i;
    % output:
    G_SE3_wam_body_{i} = G_SE3_wam_body_i;
    
    % -> compute jacobian:
    J_body_i = Lie.inv_Ad_SE3_from_SE3(G_prev) * xi_R6_0_{N_jnts-i+1};
    % output:
    J_body = [J_body_i, J_body];
end


%% [ PLOT ]:
helper.newFigure(-1);
% tiledlayout(2,1);
ax_1_1 = nexttile();

% plot summit:
utils.plot_Summit(common.SUMMIT_INIT_POSE,common.SUMMIT_INIT_POSE,'S',ax_1_1);

% plot wam base link:
utils.plot_link( ...
        G_SE3_summit_wam, ...
        eye(4), ...
        common.WAM(1).tip_frame_transformation.q', ...
        common.WAM(1).base_frame, ...
        'k', ax_1_1);

% plot rest links:
for i=1:N_jnts
    % -> plot:
    utils.plot_link( ...
        G_SE3_wam_spatial_{i}, ...
        G_SE3_0_{i}, ...
        common.WAM(i+1).tip_frame_transformation.q', ...
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
% spacial jacobian:
% Jacobian_spatial = horzcat(J_spatial_{:});
% spacial velocity:
% --- 

% --- 
helper.endSection(AUTO_CLOSE);
%% [EOF] ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 

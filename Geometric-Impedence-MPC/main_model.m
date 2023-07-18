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

% Given Angles:
JOINT_ANGLE = ones(1,N_jnts) * pi/4;
% compute FK:
exp_xi_theta_in_SE3_ = OpenChain.batch_screw_SE3_from_twist_angle_R6xR(xi_R6_0_, JOINT_ANGLE);

[G_SE3_wam_spatial_0_, J_spatial_] = OpenChain.compute_Spatial_from_SE3(G_SE3_summit_wam, xi_R6_0_, exp_xi_theta_in_SE3_);
[G_SE3_wam_body_0_, J_body_] = OpenChain.compute_Body_from_SE3(G_SE3_0_{N_jnts}, xi_R6_0_, exp_xi_theta_in_SE3_);

% [G_SE3_wam_spatial_, J_spatial_] = OpenChain.compute_Spatial(G_SE3_summit_wam, xi_R6_0_, JOINT_ANGLE);
% [G_SE3_wam_body_, J_body_] = OpenChain.compute_Body(G_SE3_0_{N_jnts}, xi_R6_0_, JOINT_ANGLE);


%% Compute Moments and Inertia:



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
% spacial jacobian:
% Jacobian_spatial = horzcat(J_spatial_{:});
% spacial velocity:
% --- 

% --- 
helper.endSection(AUTO_CLOSE);
%% [EOF] ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 

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
helper.setLogLevel("error")         % [ all, debug, error, info ]
validate.set_validatorLevel("all")  % [ all, none ]
% --- 
helper.endSection(AUTO_CLOSE);
%% Initial g(0)) ===== ===== ===== ===== ===== ===== =====:
DIR = helper.declareSection("test", "init_wam", SAVE_CONSOLE, CLEAR_OUTPUT, CLOSE_WINDOW);
% --- 

% Summit Fixed Pose: [Default: eye(4)]
SUMMIT_POSE_SE3 = eye(4) * common.SUMMIT_INIT_POSE; % summit stationary --> WAM stationary

% Summit ---> WAM base frame:
R_SO3_summit = Lie.rodrigues_SO3_from_unit_R3xR([0;0;1], 0); % -> (axis,angle) to SO3:
G_SE3_summit_wam = [R_SO3_summit, common.dP_SUMMIT_WAM.'; zeros(1,3), 1]; % relative RBT

% define base spatial frame for WAM:
WAM_Spatial_0 = SUMMIT_POSE_SE3 * G_SE3_summit_wam;
% WAM_Spatial_0 = eye(4);

% [ i-1-->i: i's Links with i's output frame]
N_links = length(common.WAM);
N_jnts = N_links; % we cosider (i-1)===>(i)'s joints/frames
valid_jnts_filter = eye(N_jnts); % <--- assume all joints active by default

% [placeholder]:
joint_limits_ = cell(1,N_links); 
G_SE3_ = cell(1,N_links); 
% spatial:
xi_R6_s_ = cell(1,N_links); 
G_SE3_s_ = cell(1,N_links);
% moment:
M_R6x6_ = cell(1,N_links); 
Ad_inv_g0_s_l_ = cell(1,N_links); 
M_R6x6_spatial_ = cell(1,N_links); 

% [iterating links forward]:
for i = 1:N_links
    helper.loginfo(sprintf("> Indexing Links @ %d-->%d",i-1,i));
    % 1. obtain relative params for current link (i-1)-->(i)
    q_R3_i = common.WAM(i).tip_frame_transformation.q';
    w_R3_i = common.WAM(i).tip_frame_transformation.w';
    t_R_i  = common.WAM(i).tip_frame_transformation.t;
    w_joint_axis_i = common.WAM(i).tip_frame_revolute_joint.axis';

    % 2. output joint properties:
    joint_name_i = common.WAM(i).tip_frame_revolute_joint.name;
    joint_limits_{i} = common.WAM(i).tip_frame_revolute_joint.limits;

    % - store joint info:
    if joint_name_i == "NotActive"
        valid_jnts_filter(i,i) = 0; % invalid angle inputs
    end

    % 3. (axis,angle) to SO3:
    R_SO3_i = Lie.rodrigues_SO3_from_unit_R3xR(w_R3_i, t_R_i);
    G_SE3_{i} = [R_SO3_i, q_R3_i; zeros(1,3), 1]; % relative RBT
    
    % 4. SE3 operation:
    % - propagate base frame TF (i-1) to tip frame (i) in current link:
    if i == 1
        G_SE3_s_{i} = WAM_Spatial_0 * G_SE3_{i}; % <--- initial wam base joint tip frame (J1)
    else
        G_SE3_s_{i} = G_SE3_s_{i-1} * G_SE3_{i}; % <--- Propagation
    end
    % - obtain global representation of the joint rotation axis
    w_R3_s_i = G_SE3_s_{i}(1:3,1:3) * w_joint_axis_i; 
    
    % - twist:
    xi_R6_s_{i} = RBT.twist_R6_from_axis_point(w_R3_s_i, G_SE3_s_{i}(1:3,4));
    
%     % [ spatial twist ] 
%     helper.logdebug(helper.a2str("xi_R6_s_{i}",xi_R6_s_{i}));
%     % [ spatial transformation ] 
%     helper.logdebug(helper.a2str("G_SE3_s_{i}",G_SE3_s_{i}));
end

xi_R6_s = cat(2,xi_R6_s_{:})

%% [INIT Jnt Angles]:
% --- 
JOINT_ANGLE = zeros(1,N_jnts)';
% JOINT_ANGLE = pi/4 * valid_jnts_filter * ones(1,N_jnts)';
%JOINT_ANGLE(3) = 0.2;

%% [iterating links backward]:
G_SE3_rvr_ = cell(1,N_links);
G_SE3_b_ = cell(1,N_links);
xi_R6_b_ = cell(1,N_links);

G_TOOL_0 = eye(4);

xi_R6_b_{N_links} = [0,0,0,0,0,1]';
% todo: please make sure the produced link can provide same jacobian
for i = N_links:-1:2
    helper.loginfo(sprintf("> Indexing Links @ %d-->%d",i+1,i));
    % 1. obtain relative params for current link (i)-->(i+1)
    G_SE3_rvr_{i} = RBT.inverse_SE3(G_SE3_{i}); % relative RBT inversed
    % - propagate base frame inverse TF (i+1) to tip frame (i) in current link:
    if i == N_links
        G_SE3_b_{i} =  G_SE3_rvr_{i} * G_TOOL_0;
    else
        G_SE3_b_{i} =  G_SE3_b_{i+1} * G_SE3_rvr_{i}; % <--- Propagation
    end

    % - obtain global representation of the joint rotation axis
    w_joint_axis_i = common.WAM(i).tip_frame_revolute_joint.axis';
    w_R3_b_i = G_SE3_b_{i}(1:3,1:3) * w_joint_axis_i; 
    
    % - twist:
    xi_R6_b_{i-1} = RBT.twist_R6_from_axis_point(w_R3_b_i, G_SE3_b_{i}(1:3,4));

    
%     % - compute generalized inertia matrix: ---->O---x--->---
%     m_R_i = common.WAM(i).mass;
%     mc_R3_i = common.WAM(i).tip_frame_mass_center';
%     I_mc_R3x3_i = common.WAM(i).tip_frame_MoI_at_mass; % inertia at mass center , aligned with output frame axis
%     mc_R3_i_hat = Lie.hat_so3_from_R3(mc_R3_i);
%     % - inertia at mass center , aligned with output frame axis:
%     % M_mass_center_i = [
%     %     m_R_i * eye(3)  , zeros(3,3);
%     %     zeros(3,3)      , I_mc_R3x3_i;
%     % ];
%     % - translating inertia matrix from center-of-mass to the output frame:
%     M_R6x6_{i} = [
%         m_R_i * eye(3)       , -m_R_i*mc_R3_i_hat; ...
%         m_R_i * mc_R3_i_hat  , I_mc_R3x3_i - m_R_i * mc_R3_i_hat^2
%     ]; % [pg 288, ]
% 
%     % - adjoint (0)
%     % Ad_inv_g0_s_l_{i} = Lie.Ad_SE3_from_SE3(RBT.inverse_SE3(G_SE3_s_i));
%     Ad_inv_g0_s_l_{i} = Lie.inv_Ad_SE3_from_SE3(G_SE3_s_i); % equivalent
% 
%     % - Inertia of the ith link reflected into the base spatial frame:
%     % [4.28, Murray]
%     M_R6x6_spatial_{i} = Ad_inv_g0_s_l_{i}' * M_R6x6_{i} * Ad_inv_g0_s_l_{i}; 
end

xi_R6_b = cat(2,xi_R6_b_{:})
% inv_ad_gst = Lie.inv_Ad_SE3_from_SE3(G_SE3_wam_spatial_ours);

%% Body) ===== ===== ===== ===== ===== ===== =====:
DIR = helper.declareSection("test", "body", SAVE_CONSOLE, CLEAR_OUTPUT, CLOSE_WINDOW);


% --- 
% Body FK:
% --- 
Blist = cat(2, xi_R6_b_{:});
% G_SE3_wam_body_MR = OpenChainMR.body_forward_kinematics(G_SE3_b_{1}, Blist, JOINT_ANGLE)
% G_SE3_wam_body_MR_2 = G_SE3_b_{1} * OpenChainMR.spatial_forward_kinematics(eye(4), Blist, JOINT_ANGLE) 
% % [NOTE]: as expected, they are equivalent
% validate.if_equivalent(G_SE3_wam_body_MR, G_SE3_wam_body_MR_2, "MR: spatial vs body function");
% % equivalent ours:
% exp_xi_theta_in_SE3_body_ = OpenChain.batch_screw_SE3_from_twist_angle_R6xR(xi_R6_b_, JOINT_ANGLE);
% G_SE3_wam_body_MR_ours = OpenChain.compute_Body_FK(exp_xi_theta_in_SE3_body_,G_SE3_b_{1})
% % Test:
% validate.if_equivalent(G_SE3_wam_body_MR, G_SE3_wam_body_MR_ours, "MR vs our function");
% % more:
% G_SE3_wam_body_MR_ours_ = OpenChain.compute_Body_FK_for_all_joints(exp_xi_theta_in_SE3_body_, G_SE3_b_);
% % Test:
% validate.if_equivalent(G_SE3_wam_body_MR_ours_{1}, G_SE3_wam_body_MR_ours, "our batch body function");

% --- 
% Body Jacobian:
% --- 
% MR:
J_wam_body_MR = OpenChainMR.body_jacobian(Blist, JOINT_ANGLE)
% ours:
neg_exp_xi_theta_in_SE3_body_ = OpenChain.batch_screw_SE3_from_twist_angle_R6xR(xi_R6_b_, - JOINT_ANGLE);
J_wam_body_ours = OpenChain.compute_Body_Jacobian(xi_R6_b_, neg_exp_xi_theta_in_SE3_body_)
% Test:
validate.if_equivalent(J_wam_body_MR, J_wam_body_ours, "compute_Body_Jacobian");


%% Spatial) ===== ===== ===== ===== ===== ===== =====:
DIR = helper.declareSection("test", "spatial", SAVE_CONSOLE, CLEAR_OUTPUT, CLOSE_WINDOW);
% --- 

% --- 
% spacial FK:
% --- 
Slist = cat(2, xi_R6_s_{:});
% mr:
G_SE3_wam_spatial_MR = OpenChainMR.spatial_forward_kinematics(G_SE3_s_{end}, Slist, JOINT_ANGLE);
% equivalent ours:
exp_xi_theta_in_SE3_ = OpenChain.batch_screw_SE3_from_twist_angle_R6xR(xi_R6_s_, JOINT_ANGLE);
G_SE3_wam_spatial_ours = OpenChain.compute_Spatial_FK(exp_xi_theta_in_SE3_, G_SE3_s_{end});
% check:
validate.if_equivalent(G_SE3_wam_spatial_MR, G_SE3_wam_spatial_ours, "Spatial Forward Kinematics");
% compute all joints FK:
G_SE3_wam_spatial_ours_ = OpenChain.compute_Spatial_FK_for_all_joints(exp_xi_theta_in_SE3_, G_SE3_s_);

% --- 
% spacial jacobian:
% --- 
% MR:
J_wam_spatial_MR = OpenChainMR.spatial_jacobian(Slist, JOINT_ANGLE);
% ours:
J_wam_spatial_ours = OpenChain.compute_Spatial_Jacobian(xi_R6_s_, exp_xi_theta_in_SE3_);
validate.if_equivalent(J_wam_spatial_MR, J_wam_spatial_ours, "Spatial Jacobian")

% --- 
% body jacobian:
% --- 
ad_gst = Lie.Ad_SE3_from_SE3(G_SE3_wam_spatial_ours);
% ours w/ spatial:
J_wam_body_ours_v2 = OpenChain.compute_Body_Jacobian_from_Spatial(xi_R6_s_, exp_xi_theta_in_SE3_, G_SE3_s_{end})
J_wam_spatial_from_body = ad_gst * J_wam_body_ours_v2; % [Murray 3.56]
validate.if_equivalent(J_wam_spatial_from_body, J_wam_spatial_ours, "J_spatial = Ad_gst * J_body [Murray 3.56]")

validate.if_equivalent(J_wam_body_ours_v2, J_wam_body_ours, "J_body_from_spatial = J_body")




%% PLOT) ===== ===== ===== ===== ===== ===== =====:
helper.endSection(AUTO_CLOSE);
DIR = helper.declareSection("test", "plot_wam", SAVE_CONSOLE, CLEAR_OUTPUT, CLOSE_WINDOW);
figure(1)
% helper.newFigure(-1);
% tiledlayout(2,1);
ax_1_1 = nexttile();

% [ Spatial ]
% plot summit:
utils.plot_Summit(SUMMIT_POSE_SE3,SUMMIT_POSE_SE3,'S',ax_1_1);

% plot base link:
utils.plot_link( ...
    WAM_Spatial_0, ...
    common.WAM(1).tip_frame_transformation.q', ...
    common.WAM(1).base_frame, ...
    'k','--');
% plot joint links:
for i=2:N_links
    utils.plot_link( ...
        G_SE3_wam_spatial_ours_{i-1}, ...
        common.WAM(i).tip_frame_transformation.q', ...
        common.WAM(i).base_frame, ...
        'k','-.');
end
% plot tip frame:
utils.plot_link( ...
    G_SE3_wam_spatial_ours_{N_links}, ...
    [0,0,0]', ...
    common.WAM(N_links).tip_frame, ...
    'k','-');

% [ Body ]
% % plot joint links:
% for i=N_links:-1:N_links-1
%     utils.plot_link( ...
%         G_SE3_wam_spatial_ours_{i-1}, ...
%         common.WAM(i).tip_frame_transformation.q', ...
%         common.WAM(i).base_frame, ...
%         'k','-.');
% end

grid on;
axis equal;
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')
% axis([-1 1 -1 1 0 1])
% helper.saveFigure([400,600], DIR, "FK")

% --- 
helper.endSection(AUTO_CLOSE);
%% [EOF] ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 

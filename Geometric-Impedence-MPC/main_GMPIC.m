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
DIR = helper.declareSection("gmpic", "init_wam", SAVE_CONSOLE, CLEAR_OUTPUT, CLOSE_WINDOW);
% --- 

% --- 
% Initialization:
% --- 
[WAM_Spatial_0, xi_R6_s_, G_SE3_s_, M_R6x6_, G_SE3_, joint_limits_, valid_jnts_filter] = model_initialization();
% convert to MR:
Slist = cat(2, xi_R6_s_{:});


% given JOINT_ANGLE:
% JOINT_ANGLE = zeros(1,N_jnts)';
N_jnts = length(valid_jnts_filter);
JOINT_ANGLE = pi/4 * valid_jnts_filter * ones(1,N_jnts)';

% --- 
% Jacobi:
% --- 
J_wam_spatial_MR = OpenChainMR.spatial_jacobian(Slist, JOINT_ANGLE);

% --- 
% mass matrix:
% --- 
Mlist = cat(3, eye(4), G_SE3_{:});
Glist = cat(3, M_R6x6_{:});
mass_MR = OpenChainMR.mass_matrix(JOINT_ANGLE, Mlist, Glist, Slist)
gravity_MR = OpenChainMR.gravity_forces(JOINT_ANGLE, [0,0,9.81], Mlist, Glist, Slist)

% --- 
% FK SIM:
% --- 
exp_xi_theta_in_SE3_ = OpenChain.batch_screw_SE3_from_twist_angle_R6xR(xi_R6_s_, JOINT_ANGLE);
G_SE3_wam_spatial_ours_ = OpenChain.compute_Spatial_FK_for_all_joints(exp_xi_theta_in_SE3_, G_SE3_s_);
helper.newFigure(-1);
plot_WAM(true, G_SE3_wam_spatial_ours_, WAM_Spatial_0);




% --- 
helper.endSection(AUTO_CLOSE);
%% [EOF] ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 
function [WAM_Spatial_0, xi_R6_s_, G_SE3_s_, M_R6x6_, G_SE3_, joint_limits_, valid_jnts_filter] = model_initialization()
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
    
        %%% compute generalized inertia matrix: ---->O---x--->---
        % - fectch mass, com, MoI:
        m_R_i = common.WAM(i).mass;
        q_mc_R3_i = common.WAM(i).tip_frame_mass_center';
        I_R3x3_i = common.WAM(i).tip_frame_MoI_at_mass; 
        
        M_R6x6_{i} = OpenChain.init_link_generalized_inertia_matrix(m_R_i, I_R3x3_i);
    end
    
    xi_R6_s = cat(2,xi_R6_s_{:})
end

function plot_WAM(if_include_summit, G_SE3_wam_spatial_, WAM_Spatial_0)
    % [ Spatial ]
    % plot summit:
    if if_include_summit
        SUMMIT_POSE_SE3 = eye(4) * common.SUMMIT_INIT_POSE;
        utils.plot_Summit(SUMMIT_POSE_SE3,SUMMIT_POSE_SE3,'S');
    end
    
    % plot base link:
    utils.plot_link( ...
        WAM_Spatial_0, ...
        common.WAM(1).tip_frame_transformation.q', ...
        common.WAM(1).base_frame, ...
        'k','--');
    % plot joint links:
    N_links = length(G_SE3_wam_spatial_);
    for i=2:N_links
        utils.plot_link( ...
            G_SE3_wam_spatial_{i-1}, ...
            common.WAM(i).tip_frame_transformation.q', ...
            common.WAM(i).base_frame, ...
            'k','-.');
    end
    % plot tip frame:
    utils.plot_link( ...
        G_SE3_wam_spatial_{N_links}, ...
        [0,0,0]', ...
        common.WAM(N_links).tip_frame, ...
        'k','-');

    grid on;
    axis equal;
    xlabel('X Axis')
    ylabel('Y Axis')
    zlabel('Z Axis')
end
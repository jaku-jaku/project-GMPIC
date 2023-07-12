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
N_jnts = length(common.RELATIVE_Axes);
for i = 1:N_jnts
    helper.loginfo(sprintf("> Indexing Links @ %d-->%d",i-1,i));
    % - obtain relative params:
    q_R3_i = common.RELATIVE_Links(i,:).';
    w_R3_i = common.RELATIVE_Axes(i,:).';
    t_R_i  = common.RELATIVE_Axes_Rotations(i).';
    
    % -> (axis,angle) to SO3:
    R_SO3_i = rodrigues_SO3_from_R3xR(w_R3_i, t_R_i);

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
    w_R3_0_i = G_SE3_0_i(1:3,1:3) * JOINT_CAD_AXIS; 
    
    % - twist:
    xi_R6_0_i = [q_R4_0_i(1:3); w_R3_0_i];
    helper.logdebug(helper.a2str("xi_R6_0_i",xi_R6_0_i));
    
    % (cache):
    xi_R6_0_{i} = xi_R6_0_i;
    G_SE3_0_{i} = G_SE3_0_i; % transformation
    helper.logdebug(helper.a2str("G_SE3_0_i",G_SE3_0_i));
end

% [ Define ]

% --- 
helper.endSection(AUTO_CLOSE);
%% Plot) ===== ===== ===== ===== ===== ===== =====:
DIR = helper.declareSection("test", "plot_wam", SAVE_CONSOLE, CLEAR_OUTPUT, CLOSE_WINDOW);
% --- 
helper.newFigure(-1);
i = 1;
% tiledlayout(2,1);
ax_1_1 = nexttile();

% -> (axis,angle) to SO3:
R_SO3_summit = rodrigues_SO3_from_R3xR([0;0;1], 0);
R_SO3_wam = rodrigues_SO3_from_R3xR([0;0;1], 0);
% (homogeneous):
G_SE3_summit = [R_SO3_summit, common.G_SUMMIT_WAM.'; zeros(1,3), 1]; % relative RBT
temp_mat_w_s = common.SUMMIT_INIT_POSE;

utils.plot_Summit(temp_mat_w_s,common.SUMMIT_INIT_POSE,'S',ax_1_1);

INIT_JOINT_HOME = [0 0 0 0 0 0 0];
for i=1:N_jnts
    theta = INIT_JOINT_HOME(i);
    if i == 1
        temp_mat_w_{1} = temp_mat_w_s * G_SE3_summit;
    else
        G_i_t_{i} = screw_SE3_from_xi_theta(xi_R6_0_{i}, theta);
        temp_mat_w_{i} = temp_mat_w_{i-1} * G_i_t_{i};
    end

    utils.plot_link( ...
        temp_mat_w_{i}, ...
        G_SE3_0_{i}, ...
        common.RELATIVE_Links(i+1,:).', ...
        sprintf('%d',i-1),'k',ax_1_1);
end
grid on;
axis equal;
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')

% axis([-1 1 -1 1 0 1])

% --- 
helper.endSection(AUTO_CLOSE);
%% [EOF] ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 

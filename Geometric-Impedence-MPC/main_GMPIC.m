close all;
clear all;
clc;

%% [INIT] ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 
% --- 
addpath(genpath("utils"));
% --- 
helper.logtitle("Initialize");
% --- 
% - USER PARAMs:
SAVE_CONSOLE = false;
CLEAR_OUTPUT = false;
CLOSE_WINDOW = true; % pre-closing
AUTO_CLOSE   = false;

helper.createFolder("output/test", false);
helper.setLogLevel("error")         % [ all, debug, error, info ]
validate.set_validatorLevel("all")  % [ all, none ]

VIEW_BNDRY = [-1 1 -1 1 0 2];
VIEW_ANGLE = [45 10];
VIEW_DIMENSION = [400,400];
% --- 
% - Simulation Params:
FB_Ctrl_CONFIG = struct("Kp", 0.8, "Ki", 0.8, "Kd", 1);
MODES_EE_F = ["Zero_F_EE","Const_F_EE","Dynamic_F_EE"];
MODES_TRAJ = ["Steady","Const_dTheta","Dynamic_dTheta"];

MODE = [MODES_EE_F(3), MODES_TRAJ(3)];
Exp_Title = strcat(MODE(:));

% - compute tag:
FB_CTRL_TAG = sprintf("P[%s]I[%s]D[%s]",...
    helper.float2str(FB_Ctrl_CONFIG.Kp),...
    helper.float2str(FB_Ctrl_CONFIG.Ki),...
    helper.float2str(FB_Ctrl_CONFIG.Kd));
% --- 
helper.endSection(AUTO_CLOSE);
%% Initial g(0)) ===== ===== ===== ===== ===== ===== =====:
DIR = helper.declareSection("gmpic", Exp_Title, SAVE_CONSOLE, CLEAR_OUTPUT, CLOSE_WINDOW);
% --- 
helper.logtitle("WAM INIT")
helper.newFigure(-1);
c0 = WAM_model.model_init_zero_config(); 
WAM_model.plot_zero_config(c0);
utils.plot3DBoundary(VIEW_BNDRY);
view(VIEW_ANGLE)
axis equal;
helper.saveFigure(VIEW_DIMENSION, DIR, "model");


%% gen angles:
helper.logtitle("Generate Joint Trajectory Ref")
N_end = 101; % total steps
list_of_thetas = zeros(8,N_end);
list_of_force_at_EE = zeros(6,N_end);
FPS = 10;
dT = 1/FPS; %s
for i=2:N_end 
    switch MODE(2)
        case "Steady"
            list_of_thetas(:,i) = [pi/2, -pi/2, 0, 0, 0, 0, 0, 0 ]';
        case "Const_dTheta"
            list_of_thetas(:,i) = [-2.6+pi/50*i, -pi/100*i, 0, 0, 0, 0, 0, 0 ]';
        case "Dynamic_dTheta"
            list_of_thetas(:,i) = [-2.6+2*pi*sin(i/100*2*pi), -pi/100*i, 0, 0, 0, 0, 0, 0 ]';
    end
    % apply input constraints:
    list_of_thetas(:,i) = WAM_model.filter_angles(c0,list_of_thetas(:,i));
    list_of_thetas(:,i) = WAM_model.joint_constraints(c0,list_of_thetas(:,i));
    % inject disturbance at the EE:
    switch MODE(1)
        case "Zero_F_EE"
            continue % do nothing
        case "Const_F_EE"
            list_of_force_at_EE(:,i) = ones(6,1);
        case "Dynamic_F_EE"
            list_of_force_at_EE(:,i) = ones(6,1)*sin(pi*i/10);
    end
end
list_of_thetas(:,1)=list_of_thetas(:,2);
%% compute fwd traj:
helper.logtitle("Compute FWD Traj")
spatial_frame_SE3_ = WAM_model.compute_EE_trajectory_with(c0,list_of_thetas);

% plot:
helper.newFigure(-1);
utils.plot_trajectory_from_SE3(spatial_frame_SE3_,0.1, "-", false, 10);
utils.plot3DBoundary(VIEW_BNDRY);
WAM_model.plot_at(c0, list_of_thetas(:,1));
view(VIEW_ANGLE)
axis equal;
helper.saveFigure(VIEW_DIMENSION, DIR, "reference_trajectory");


%% compute trajectory after feedback linearizd
% ation controller with dynamics:
helper.logtitle("Compute Dynamics Controller")
data = WAM_model.dynamic_sim(c0, list_of_thetas(:,1), list_of_thetas(:,2:end), ...
    list_of_force_at_EE, dT, FB_Ctrl_CONFIG);
% compute fwd traj:
spatial_frame_SE3_fb_ = WAM_model.compute_EE_trajectory_with(c0, data.theta);
% plot:
helper.newFigure(-1);
utils.plot_trajectory_from_SE3(spatial_frame_SE3_,0.1, "--", false, 10);
utils.plot_trajectory_from_SE3(spatial_frame_SE3_fb_,0.1, "-", false, 10);
utils.plot3DBoundary(VIEW_BNDRY);
WAM_model.plot_at(c0, list_of_thetas(:,1));
view(VIEW_ANGLE)
axis equal;
file_name = strcat("trajectory_fb_controller_", FB_CTRL_TAG);
helper.saveFigure(VIEW_DIMENSION, DIR, file_name);

%% plot joints:
helper.newFigure(-1);
utils.plot_position_vs_t(...
    {spatial_frame_SE3_,spatial_frame_SE3_fb_}, ...
    ["ref", "PID"], ...
    dT);
helper.saveFigure(VIEW_DIMENSION, DIR, strcat("EE_position_",FB_CTRL_TAG));

%% animate:
pause(0.00001)
helper.logtitle("Animating ...")
WAM_model.animate_with(c0, data.theta, VIEW_DIMENSION, VIEW_BNDRY, VIEW_ANGLE, DIR, FB_CTRL_TAG, FPS);

% --- 
helper.endSection(AUTO_CLOSE);
%% [EOF] ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 

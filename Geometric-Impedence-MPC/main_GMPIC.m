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
Exp_Title = "Const_F_EE";
% Exp_Title = "Zero_F_EE";

% --- 
helper.endSection(AUTO_CLOSE);
%% Initial g(0)) ===== ===== ===== ===== ===== ===== =====:
DIR = helper.declareSection("gmpic", Exp_Title, SAVE_CONSOLE, CLEAR_OUTPUT, CLOSE_WINDOW);
% --- 
helper.logtitle("WAM INIT")
model = WAM_model();
helper.newFigure(-1);
model.plot_now();
utils.plot3DBoundary(VIEW_BNDRY);
view(VIEW_ANGLE)
axis equal;
helper.saveFigure(VIEW_DIMENSION, DIR, "model");

% model = model.set_angles([ 0, -2, 0, 3.13, 0, 0, 0 , 0]');
% model.plot_now();

%% gen angles:
helper.logtitle("Generate Joint Trajectory Ref")
N_end = 101; % total steps
list_of_thetas = zeros(8,N_end);
list_of_force_at_EE = zeros(6,N_end);
FPS = 10;
dT = 1/FPS; %s
for i=2:N_end 
    list_of_thetas(:,i) = [-2.6+pi/50*i, -pi/100*i, 0, 0, 0, 0, 0, 0 ]';
    % apply input constraints:
    list_of_thetas(:,i) = model.filter_angles(list_of_thetas(:,i));
    list_of_thetas(:,i) = model.joint_constraints(list_of_thetas(:,i));
    % inject disturbance at the EE:
    if Exp_Title == "Zero_F_EE"
        continue % do nothing
    elseif Exp_Title == "Const_F_EE"
        list_of_force_at_EE(:,i) = ones(6,1);
    elseif Exp_Title == "Dynamic_F_EE"
        list_of_force_at_EE(:,i) = ones(6,1)*sin(pi*i/10);
    end
end
list_of_thetas(:,1)=list_of_thetas(:,2);
model.set_angle(list_of_thetas(:,1));
%% compute fwd traj:
helper.logtitle("Compute FWD Traj")
spatial_frame_SE3_ = model.compute_EE_trajectory_with(list_of_thetas);

% plot:
helper.newFigure(-1);
utils.plot_trajectory_from_SE3(spatial_frame_SE3_,0.1, "-", false, 10);
utils.plot3DBoundary(VIEW_BNDRY);
model.plot_now();
view(VIEW_ANGLE)
axis equal;
helper.saveFigure(VIEW_DIMENSION, DIR, "trajectory");


%% compute trajectory after feedback linearization controller with dynamics:
helper.logtitle("Compute Dynamics Controller")
FB_Ctrl_CONFIG = struct("Kp", 0.8, "Ki", 0.8, "Kd", 1);
data = model.dynamic_sim(list_of_thetas(:,1), list_of_thetas(:,2:end), ...
    list_of_force_at_EE, dT, FB_Ctrl_CONFIG);
% compute fwd traj:
spatial_frame_SE3_fb_ = model.compute_EE_trajectory_with(data.theta);
% plot:
helper.newFigure(-1);
utils.plot_trajectory_from_SE3(spatial_frame_SE3_,0.1, "--", false, 10);
utils.plot_trajectory_from_SE3(spatial_frame_SE3_fb_,0.1, "-", false, 10);
utils.plot3DBoundary(VIEW_BNDRY);
model.plot_now();
view(VIEW_ANGLE)
axis equal;
TAG = sprintf("P[%s]I[%s]D[%s]",...
    helper.float2str(FB_Ctrl_CONFIG.Kp),...
    helper.float2str(FB_Ctrl_CONFIG.Ki),...
    helper.float2str(FB_Ctrl_CONFIG.Kd));
file_name = strcat("trajectory_fb_controller_",TAG);
helper.saveFigure(VIEW_DIMENSION, DIR, file_name);

%% plot joints:
helper.newFigure(-1);
utils.plot_position_vs_t(...
    {spatial_frame_SE3_,spatial_frame_SE3_fb_}, ...
    ["ref", "PID"], ...
    dT);
helper.saveFigure(VIEW_DIMENSION, DIR, strcat("EE_position",TAG));

%% animate:
pause(0.00001)
helper.logtitle("Animating ...")
model.animate_with(data.theta, VIEW_DIMENSION, VIEW_BNDRY, VIEW_ANGLE, DIR, TAG, FPS);

% --- 
helper.endSection(AUTO_CLOSE);
%% [EOF] ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 

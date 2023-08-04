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
% --- 
helper.endSection(AUTO_CLOSE);
%% Initial g(0)) ===== ===== ===== ===== ===== ===== =====:
DIR = helper.declareSection("gmpic", "init_wam", SAVE_CONSOLE, CLEAR_OUTPUT, CLOSE_WINDOW);
% --- 

model = WAM_model();
helper.newFigure(-1);
model.plot_now();
utils.plot3DBoundary(VIEW_BNDRY);
view(VIEW_ANGLE)
axis equal;
helper.saveFigure([600,300], DIR, "model");

% model = model.set_angles([ 0, -2, 0, 3.13, 0, 0, 0 , 0]');
% model.plot_now();

%% gen angles:
list_of_thetas = zeros(8,100);
for i=1:100
    list_of_thetas(:,i) = [ pi/100*i, -pi/200*i, pi/100*i, -pi/100*i, pi/100*i, 0, 0, 0 ]';
end
%% plot traj:
helper.newFigure(-1);
spatial_frame_SE3_ = model.compute_EE_trajectory_with(list_of_thetas);
utils.plot_trajectory_from_SE3(spatial_frame_SE3_, 0.1, "-", false);
utils.plot3DBoundary(VIEW_BNDRY);
model.plot_now();
view(VIEW_ANGLE)
axis equal;
helper.saveFigure([600,300], DIR, "trajectory");

%% animate:
model.animate_with(list_of_thetas, VIEW_BNDRY, VIEW_ANGLE, DIR, "test", 10)

% --- 
helper.endSection(AUTO_CLOSE);
%% [EOF] ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 

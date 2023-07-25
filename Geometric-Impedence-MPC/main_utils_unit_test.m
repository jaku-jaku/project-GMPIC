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
CLOSE_WINDOW = true; % pre closing
AUTO_CLOSE = true;

helper.createFolder("output/test", false);

helper.setLogLevel("all")
% helper.setLogLevel("error") % Hide debug msg!!
validate.set_validatorLevel("all")
% --- 
helper.endSection(AUTO_CLOSE);
%% Tests) ===== ===== ===== ===== ===== ===== =====:
DIR = helper.declareSection("test", "1", SAVE_CONSOLE, CLEAR_OUTPUT, CLOSE_WINDOW);
% --- 
helper.logtitle("hat, vee")
vec_vw = [4; 5; 6; 1; 2; 3];
mat_se3_target =[
    0    -3     2     4;
    3     0    -1     5;
   -2     1     0     6;
    0     0     0     0
];
mat_se3 = Lie.hat_se3_from_R6(vec_vw)
validate.test_compare(mat_se3, mat_se3_target, "hat_se3_from_R6");

vec_se3 = Lie.vee_R6_from_se3(mat_se3_target)'
validate.test_compare(vec_se3, vec_vw, "vee_R6_from_se3");



helper.logtitle("rodrigue")
se3mat = [ 0,      0,       0,      0;
           0,      0, -1.5708, 2.3562;
           0, 1.5708,       0, 2.3562;
           0,      0,       0,      0];
T_target = [
   1.0000         0         0         0;
        0    0.0000   -1.0000   -0.0000;
        0    1.0000    0.0000    3.0000;
        0         0         0    1.0000;
];

T = Lie.exp_map_SE3_from_R6(Lie.vee_R6_from_se3(se3mat)')
validate.test_compare(T, T_target, "exp_map_SE3_from_R6");

T = [[1, 0, 0, 0]; [0, 0, -1, 0]; [0, 1, 0, 3]; [0, 0, 0, 1]];
T_target = [
 0         0         0         0;
 0         0   -1.5708    2.3562;
 0    1.5708         0    2.3562;
 0         0         0         0;
];
T_mat = Lie.log_se3_from_SE3(T)
validate.test_compare(T_mat, T_target, "log_se3_from_SE3");

%% ---
helper.logtitle("spatial_forward_kinematics")

M = [[-1, 0, 0, 0]; [0, 1, 0, 6]; [0, 0, -1, 2]; [0, 0, 0, 1]];
Slist = [[ 4; 0;    0; 0; 0;  1], ...
         [ 0; 1;    0; 0; 0;  0], ...
         [-6; 0; -0.1; 0; 0; -1]];
thetalist =[pi / 2; 3; pi];

T_target = [
-0.0000    1.0000         0   -5.0000;
 1.0000    0.0000         0    4.0000;
      0         0   -1.0000    1.6858;
      0         0         0    1.0000;
];

T_mat = OpenChain.spatial_forward_kinematics(M, Slist, thetalist)

validate.test_compare(T_mat, T_target, "spatial_forward_kinematics");

%% ---
helper.logtitle("spatial_inverse_kinematics")

Slist = [[  4; 0;    0; 0; 0;  1], ...
         [  0; 1;    0; 0; 0;  0], ...
         [ -6; 0; -0.1; 0; 0; -1]];
M = [[-1, 0, 0, 0]; [0, 1, 0, 6]; [0, 0, -1, 2]; [0, 0, 0, 1]];
T = [[0, 1, 0, -5]; [1, 0, 0, 4]; [0, 0, -1, 1.6858]; [0, 0, 0, 1]];
thetalist0 = [1.5; 2.5; 3];
eomg = 0.01;
ev = 0.001;

[thetalist, success] = OpenChain.spatial_inverse_kinematics(...
    Slist, M, T, thetalist0, eomg, ev)

theta_target = [1.5707; 2.9997; 3.1415];
validate.test_compare(thetalist, theta_target, "spatial_inverse_kinematics");

%% ---
helper.logtitle("body_forward_kinematics")

Blist = [[ 2; 0; 0; 0; 0; -1], ...
         [0; 1; 0; 0; 0; 0], ...
         [0; 0; 0.1; 0; 0; 1]];
M = [[-1, 0, 0, 0]; [0, 1, 0, 6]; [0, 0, -1, 2]; [0, 0, 0, 1]];
thetalist =[pi / 2; 3; pi];

T_target = [
-0.0000    1.0000         0   -5.0000;
 1.0000    0.0000         0    4.0000;
      0         0   -1.0000    1.6858;
      0         0         0    1.0000;
];

T_mat = OpenChain.body_forward_kinematics(M, Blist, thetalist)

validate.test_compare(T_mat, T_target, "body_forward_kinematics");
%% ---
helper.logtitle("body_inverse_kinematics")

Blist = [[ 2; 0; 0; 0; 0; -1], ...
         [0; 1; 0; 0; 0; 0], ...
         [0; 0; 0.1; 0; 0; 1]];
M = [[-1, 0, 0, 0]; [0, 1, 0, 6]; [0, 0, -1, 2]; [0, 0, 0, 1]];
T = [[0, 1, 0, -5]; [1, 0, 0, 4]; [0, 0, -1, 1.6858]; [0, 0, 0, 1]];
thetalist0 = [1.5; 2.5; 3];
eomg = 0.01;
ev = 0.001;

[thetalist, success] = OpenChain.body_inverse_kinematics(...
    Blist, M, T, thetalist0, eomg, ev)

theta_target = [1.5707; 2.9997; 3.1415];
validate.test_compare(thetalist, theta_target, "body_inverse_kinematics");


%% [OpenChain Lib Testing]:
% testing migrations from modern robotics:
Slist = [[  0; 0.2; 0.2; 0; 0; 1], ...
         [  2;   0;   3; 1; 0; 0], ...
         [  0;   2;   1; 0; 1; 0], ...
         [0.2; 0.3; 0.4; 1; 0; 0]];
theta = [0.2; 1.1; 0.1; 1.2];

% test:
helper.logtitle("Jacobian")
Js_target =    [
         0    1.9522   -2.2164   -0.5116;
    0.2000    0.4365   -2.4371    2.7754;
    0.2000    2.9603    3.2357    2.2251;
         0    0.9801   -0.0901    0.9575;
         0    0.1987    0.4446    0.2849;
    1.0000         0    0.8912   -0.0453;
];
Js = OpenChain.spatial_jacobian(Slist, theta)
validate.test_compare(Js, Js_target, "Space Jacobian");

Jb_target =    [
      2.3259    1.6681    0.5641    0.2000;
     -1.4432    2.9456    1.4331    0.3000;
     -2.0664    1.8288   -1.5887    0.4000;
     -0.0453    0.9950         0    1.0000;
      0.7436    0.0930    0.3624         0;
     -0.6671    0.0362   -0.9320         0;
];
Jb = OpenChain.body_jacobian(Slist, theta)
validate.test_compare(Jb, Jb_target, "Body Jacobian");


%% Example Input (3 Link Robot):
% validate.set_validatorLevel('all')
theta = [0.1; 0.1; 0.1];
d_theta = [0.1; 0.2; 0.3];
dd_theta = [2; 1.5; 1];
g = [0; 0; -9.8];
Ftip = [1; 1; 1; 1; 1; 1];
M01 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.089159]; [0, 0, 0, 1]];
M12 = [[0, 0, 1, 0.28]; [0, 1, 0, 0.13585]; [-1, 0 ,0, 0]; [0, 0, 0, 1]];
M23 = [[1, 0, 0, 0]; [0, 1, 0, -0.1197]; [0, 0, 1, 0.395]; [0, 0, 0, 1]];
M34 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
G1 = diag([3.7, 3.7, 3.7, 0.010267, 0.010267, 0.00666]);
G2 = diag([8.393, 8.393, 8.393, 0.22689, 0.22689, 0.0151074]);
G3 = diag([2.275, 2.275, 2.275, 0.0494433, 0.0494433, 0.004095]);
Glist = cat(3, G1, G2, G3);
Mlist = cat(3, M01, M12, M23, M34); 
Slist = [[      0; 1;     0; 1; 0; 1], ...
         [ -0.089; 0;     0; 0; 1; 0], ...
         [ -0.089; 0; 0.425; 0; 1; 0]];

% test:
helper.logtitle("Dynamic Equations")
tau_target = [ 74.6962; -33.0675; -3.2306 ];
tau = OpenChain.inverse_dynamics(theta, d_theta, dd_theta, g, Ftip, Mlist, Glist, Slist)
validate.test_compare(tau, tau_target, "Inverse Dynamics")

C_target = [0.2645; -0.0551; -0.0069];
C = OpenChain.vel_qualdratic_force(theta, d_theta, Mlist, Glist, Slist)
validate.test_compare(C, C_target, "Vel Quadratic Force")

M_target = [
    22.5433   -0.3071   -0.0072;
    -0.3071    1.9685    0.4322;
    -0.0072    0.4322    0.1916; 
];
M = OpenChain.mass_matrix(theta, Mlist, Glist, Slist)
validate.test_compare(M, M_target, "Mass Matrix")

%%
% test:
helper.logtitle("Computed Torque")
theta = [0.1; 0.1; 0.1];
d_theta = [0.1; 0.2; 0.3];
e_int = [0.2; 0.2; 0.2];
g = [0; 0; -9.8];
M01 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.089159]; [0, 0, 0, 1]];
M12 = [[0, 0, 1, 0.28]; [0, 1, 0, 0.13585]; [-1, 0 ,0, 0]; [0, 0, 0, 1]];
M23 = [[1, 0, 0, 0]; [0, 1, 0, -0.1197]; [0, 0, 1, 0.395]; [0, 0, 0, 1]];
M34 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
G1 = diag([3.7, 3.7, 3.7, 0.010267, 0.010267, 0.00666]);
G2 = diag([8.393, 8.393, 8.393, 0.22689, 0.22689, 0.0151074]);
G3 = diag([2.275, 2.275, 2.275, 0.0494433, 0.0494433, 0.004095]);
Glist = cat(3, G1, G2, G3);
Mlist = cat(3, M01, M12, M23, M34); 
Slist = [[      0; 1;     0;1; 0; 1], ...
         [ -0.089; 0;     0;0; 1; 0], ...
         [ -0.089; 0; 0.425;0; 1; 0]];
r_theta = [1; 1; 1];
r_d_theta = [2; 1.2; 2];
r_dd_theta = [0.1; 0.1; 0.1];
Kp = 1.3;
Ki = 1.2;
Kd = 1.1;

% test:
tau_target = [133.0053; -29.9422; -3.0328];
tau = OpenChain.force_computed_at_EE(theta, d_theta, e_int, g, ... 
                            Mlist, Glist, Slist, ...
                            r_theta, r_d_theta, r_dd_theta, ...
                            Kp, Ki, Kd)
validate.test_compare(tau, tau_target, "Computed Torque")




%% --- 
helper.endSection(AUTO_CLOSE);
%% [EOF] ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 

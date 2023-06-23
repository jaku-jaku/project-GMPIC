close all;
clear all;
clc;

%% init helper:
addpath(genpath("utils"));
helper.createFolder("output/test", false);

%% so3 coord <-> so3
syms x1 x2 x3 
x = [x1; x2; x3];
xs = kron([1 2 3 4], x);

vs_h = hat_so3(x)
vs_hs = hat_so3(xs)
xs_h_v = vee_so3(vs_h)
xs_hs_v = vee_so3(vs_hs)

%% se3 coord <-> se3
syms v1 v2 v3 w1 w2 w3 
xi = [v1 v2 v3 w1 w2 w3].';
xis = kron([1 2 3 4], xi);

xis_h = hat_se3(xi)
xis_h = hat_se3(xis)

%% adjoint se3 coord -> 6x6xk
ad_se3_(xi)
ad_se3_(xis)

%% SE3
R = sym('R', [3 3]);
p = sym('p', [3 1]);

Rs = reshape(kron([1 2 3 4], R), 3, 3, 4);
ps = reshape(kron([1 2 3 4], p), 3,1,4);

%%
Ad_SE3(R, p)
Ad_SE3(Rs, ps)



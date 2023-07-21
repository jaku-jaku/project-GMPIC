classdef OpenChain
    methods(Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Ours & Efficient:
        function [G_SE3_wam_spatial_0_, J_spatial_] = compute_Spatial(G_SE3_0_spatial, xi_R6_0_, theta_)
            J_spatial_ = [];
            G_SE3_wam_spatial_0_ = {};
            
            % init:
            N_jnts = length(theta_);
            G_SE3_wam_spatial_0_{1} = G_SE3_0_spatial; % wam base location
            
            % Compute Spatial (Forward):
            for i=1:N_jnts
                % -> compute transformation from twist coordinates and angles applied:
                % |---- (more general):
                G_i_t_{i} = RBT.screw_SE3_from_xi_theta(xi_R6_0_{i}, theta_(i));
                % |---- (more efficient):
                % G_i_t_{i} = RBT.screw_SE3_from_axis_point_theta(w_R3_0_{i}, q_R3_0_{i}, theta_(i));
                
                % -> cumulate transformation:
                G_SE3_wam_spatial_0_{i+1} = G_SE3_wam_spatial_0_{i} * G_i_t_{i};
                
                % -> compute jacobian:
                J_spatial_i = Lie.Ad_SE3_from_SE3(G_SE3_wam_spatial_0_{i}) * xi_R6_0_{i};
                % output:
                J_spatial_ = [J_spatial_, J_spatial_i];
            end
        end
        function [G_SE3_wam_body_0_, J_body_] = compute_Body(G_SE3_tool_frame, xi_R6_0_, theta_)
            J_body_ = [];
            G_SE3_wam_body_0_ = {};
            
            % init:
            N_jnts = length(theta_);
            G_SE3_wam_body_0_{1} = G_SE3_tool_frame;
            
            %% Compute Body (Backward):
            for i=1:N_jnts
                % - compute backwards from tool to base:
                G_b_t_{i} = RBT.screw_SE3_from_xi_theta(xi_R6_0_{N_jnts-i+1}, theta_(N_jnts-i+1));
                % -> cumulate transformation:
                G_SE3_wam_body_0_{i+1} = G_SE3_wam_body_0_{i} * G_b_t_{i};
                
                % -> compute jacobian:
                J_body_i = Lie.inv_Ad_SE3_from_SE3(G_SE3_wam_body_0_{i}) * xi_R6_0_{N_jnts-i+1};
                % output:
                J_body_ = [J_body_i, J_body_];
            end
        end
        function [G_SE3_wam_spatial_0_, J_spatial_] = compute_Spatial_from_SE3(G_SE3_0_spatial, xi_R6_0_, exp_xi_theta_in_SE3_)
            J_spatial_ = [];
            G_SE3_wam_spatial_0_ = {};
            
            % init:
            N_jnts = length(exp_xi_theta_in_SE3_);
            G_SE3_wam_spatial_0_{1} = G_SE3_0_spatial; % wam base location
            
            % Compute Spatial (Forward):
            for i=1:N_jnts
                % -> cumulate transformation:
                G_SE3_wam_spatial_0_{i+1} = G_SE3_wam_spatial_0_{i} * exp_xi_theta_in_SE3_{i};
                % -> compute jacobian:
                J_spatial_i = Lie.Ad_SE3_from_SE3(G_SE3_wam_spatial_0_{i}) * xi_R6_0_{i};
                % output:
                J_spatial_ = [J_spatial_, J_spatial_i];
            end
        end
        function [G_SE3_wam_body_0_, J_body_] = compute_Body_from_SE3(G_SE3_tool_frame, xi_R6_0_, exp_xi_theta_in_SE3_)
            J_body_ = [];
            G_SE3_wam_body_0_ = {};
            
            % init:
            N_jnts = length(exp_xi_theta_in_SE3_);
            G_SE3_wam_body_0_{1} = G_SE3_tool_frame;
            
            %% Compute Body (Backward):
            for i=1:N_jnts
                % -> cumulate transformation:
                G_SE3_wam_body_0_{i+1} = G_SE3_wam_body_0_{i} * exp_xi_theta_in_SE3_{N_jnts-i+1};
                % -> compute jacobian:
                J_body_i = Lie.inv_Ad_SE3_from_SE3(G_SE3_wam_body_0_{i}) * xi_R6_0_{N_jnts-i+1};
                % output:
                J_body_ = [J_body_i, J_body_];
            end
        end
        function exp_xi_theta_in_SE3_ = batch_screw_SE3_from_twist_angle_R6xR(xi_R6_0_, theta_)
            exp_xi_theta_in_SE3_ = {};
            N_jnts = length(theta_);
            for i=1:N_jnts
                % -> compute transformation from twist coordinates and angles applied:
                % |---- (more general):
                exp_xi_theta_in_SE3_{i} = RBT.screw_SE3_from_xi_theta(xi_R6_0_{i}, theta_(i));
                % |---- (more efficient):
                % G_i_t_{i} = RBT.screw_SE3_from_axis_point_theta(w_R3_0_{i}, q_R3_0_{i}, theta_(i));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Modern Robotics:
        function G = body_FK(M, Blist, thetalist)
            % *** CHAPTER 9: TRAJECTORY GENERATION ***
            % Takes M: The home configuration (position and orientation) of the end-
            %          effector,
            %       Blist: The joint screw axes in the end-effector frame when the
            %              manipulator is at the home position,
            %       thetalist: A list of joint coordinates.
            % Returns G: The end-effector configuration resulting from forward kinematics
            %            (position and orientation) represented in SE(3).

            G = M;
            for i = 1: size(thetalist, 1)
                G = G * Lie.exp_map_SE3_from_R6(Blist(:, i) * thetalist(i));
            end
        end
        function G = spatial_FK(M, Slist, thetalist)
            G = M;
            for i = 1: size(thetalist, 1)
                G = Lie.exp_map_SE3_from_R6(Slist(:, i) * thetalist(i)) * G;
            end
        end
        function Js = space_jacobian(Slist, theta)
            % *** CHAPTER 5: VELOCITY KINEMATICS AND STATICS ***
            % Takes Slist: The joint screw axes in the space frame when the manipulator
            %              is at the home position, in the format of a matrix with the
            %              screw axes as the columns,
            %       theta: A list of joint coordinates. 
            % Returns the corresponding space Jacobian (6xn real numbers).
            % Example Input:
            % 
            % clear; clc;
            % Slist = [[  0; 0.2; 0.2; 0; 0; 1], ...
            %          [  2;   0;   3; 1; 0; 0], ...
            %          [  0;   2;   1; 0; 1; 0], ...
            %          [0.2; 0.3; 0.4; 1; 0; 0]];
            % theta = [0.2; 1.1; 0.1; 1.2];
            % Js = OpenChain.space_jacobian(Slist, theta)
            % 
            % Output:
            % Js =
            %             0    1.9522   -2.2164   -0.5116
            %        0.2000    0.4365   -2.4371    2.7754
            %        0.2000    2.9603    3.2357    2.2251
            %             0    0.9801   -0.0901    0.9575
            %             0    0.1987    0.4446    0.2849
            %        1.0000         0    0.8912   -0.0453
            
            Js = Slist;
            T = eye(4);
            for i = 2: length(theta)
                % [A generic one]:
                T = T * Lie.exp_map_SE3_from_R6(Slist(:, i - 1) * theta(i - 1));
                % [Assuming, Ai(:,i) contains a unit axis vector: w_R3]:
                % T = T * Lie.exp_map_SE3_from_R6xR(Slist(:, i - 1), theta(i - 1));
                Js(:, i) = Lie.Ad_SE3_from_SE3(T) * Slist(:, i);
            end
        end
        function Jb = body_jacobian(Blist, theta)
            % *** CHAPTER 5: VELOCITY KINEMATICS AND STATICS ***
            % Takes Blist: The joint screw axes in the end-effector frame when the
            %              manipulator is at the home position, in the format of a 
            %              matrix with the screw axes as the columns,
            %       theta: A list of joint coordinates.
            % Returns the corresponding body Jacobian (6xn real numbers).
            % Example Input:
            % 
            % Output:
            % Jb =
            %    2.3259    1.6681    0.5641    0.2000
            %   -1.4432    2.9456    1.4331    0.3000
            %   -2.0664    1.8288   -1.5887    0.4000
            %   -0.0453    0.9950         0    1.0000
            %    0.7436    0.0930    0.3624         0
            %   -0.6671    0.0362   -0.9320         0

            Jb = Blist;
            T = eye(4);
            for i = length(theta) - 1: -1: 1   
                % [A generic one]:
                T = T * Lie.exp_map_SE3_from_R6(-1 * Blist(:, i + 1) * theta(i + 1));
                % [Assuming, Ai(:,i) contains a unit axis vector: w_R3]:
                % T = T * Lie.exp_map_SE3_from_R6xR(-1 * Blist(:, i + 1), theta(i + 1));
                Jb(:, i) = Lie.Ad_SE3_from_SE3(T) * Blist(:, i);
            end
        end
        function taulist = inverse_dynamics(theta, d_theta, dd_theta, g, Ftip, Mlist, Glist, Slist)
            % Takes theta: n-vector of joint variables,
            %       d_theta: n-vector of joint rates,
            %       dd_theta: n-vector of joint accelerations,
            %       g: Gravity vector g,
            %       Ftip: Spatial force applied by the end-effector expressed in frame 
            %             {n+1},
            %       Mlist: List of link frames {i} relative to {i-1} at the home 
            %              position,
            %       Glist: Spatial inertia matrices Gi of the links,
            %       Slist: Screw axes Si of the joints in a space frame, in the format
            %              of a matrix with the screw axes as the columns.
            % Returns taulist: The n-vector of required joint forces/torques.
            % This function uses forward-backward Newton-Euler iterations to solve the 
            % equation:
            % taulist = Mlist(thetalist) * ddthetalist + c(thetalist, dthetalist) ...
            %           + g(thetalist) + Jtr(thetalist) * Ftip

            n = size(theta, 1);
            Mi = eye(4);
            Ai = zeros(6, n);
            AdTi = zeros(6, 6, n + 1);
            Vi = zeros(6, n + 1);
            Vdi = zeros(6, n + 1);
            Vdi(1:3, 1) = -g;
            AdTi(:, :, n + 1) = Lie.Ad_SE3_from_SE3(RBT.inverse_SE3(Mlist(:, :, n + 1)));
            Fi = Ftip;
            taulist = zeros(n, 1);
            for i=1: n    
                Mi = Mi * Mlist(:, :, i);
                Ai(:, i) = Lie.Ad_SE3_from_SE3(RBT.inverse_SE3(Mi)) * Slist(:, i);    
                % [A generic one]:
                AdTi(:, :, i) = Lie.Ad_SE3_from_SE3( ...
                    Lie.exp_map_SE3_from_R6(Ai(:, i) * -theta(i)) * RBT.inverse_SE3(Mlist(:, :, i)) ...
                );
                % [Assuming, Ai(:,i) contains a unit axis vector: w_R3]:
                % AdTi(:, :, i) = Lie.Ad_SE3_from_SE3( ...
                %     Lie.exp_map_SE3_from_R6xR(Ai(:, i), -theta(i)) * RBT.inverse_SE3(Mlist(:, :, i)) ...
                % );
                Vi(:, i + 1) = AdTi(:, :, i) * Vi(:, i) + Ai(:, i) * d_theta(i);
                Vdi(:, i + 1) = AdTi(:, :, i) * Vdi(:, i) ...
                                + Ai(:, i) * dd_theta(i) ...
                                + Lie.ad_se3_from_R6(Vi(:, i + 1)) * Ai(:, i) * d_theta(i);    
            end

            for i = n: -1: 1
                Fi = AdTi(:, :, i + 1)' * Fi + Glist(:, :, i) * Vdi(:, i + 1) ...
                    - Lie.ad_se3_from_R6(Vi(:, i + 1))' * (Glist(:, :, i) * Vi(:, i + 1));
                taulist(i) = Fi' * Ai(:, i);
            end
        end
        function c = vel_qualdratic_force(theta, d_theta, Mlist, Glist, Slist)
            % *** CHAPTER 8: DYNAMICS OF OPEN CHAINS ***
            % Takes theta: A list of joint variables,
            %       d_theta: A list of joint rates,
            %       Mlist: List of link frames i relative to i-1 at the home position,
            %       Glist: Spatial inertia matrices Gi of the links,
            %       Slist: Screw axes Si of the joints in a space frame, in the format
            %              of a matrix with the screw axes as the columns,
            % Returns c: The vector c(thetalist,dthetalist) of Coriolis and centripetal
            %            terms for a given thetalist and dthetalist.
            % This function calls InverseDynamics with g = 0, Ftip = 0, and 
            % ddthetalist = 0.
            % Example Input (3 Link Robot):
            %
            % Output:
            % c =
            %    0.2645
            %   -0.0551
            %   -0.0069
            
            c = OpenChain.inverse_dynamics(theta, d_theta, ...
                                zeros(size(theta, 1), 1), [0; 0; 0], ...
                                [0; 0; 0; 0; 0; 0], Mlist, Glist, Slist);
        end
        function M = mass_matrix(theta, Mlist, Glist, Slist)
            % *** CHAPTER 8: DYNAMICS OF OPEN CHAINS ***
            % Takes theta: A list of joint variables,
            %       Mlist: List of link frames i relative to i-1 at the home position,
            %       Glist: Spatial inertia matrices Gi of the links,
            %       Slist: Screw axes Si of the joints in a space frame, in the format
            %              of a matrix with the screw axes as the columns.
            % Returns M: The numerical inertia matrix M(thetalist) of an n-joint serial
            %            chain at the given configuration thetalist.
            % This function calls InverseDynamics n times, each time passing a 
            % ddthetalist vector with a single element equal to one and all other 
            % inputs set to zero. Each call of InverseDynamics generates a single 
            % column, and these columns are assembled to create the inertia matrix.
            % Example Input (3 Link Robot):
            % 
            % Output:
            % M =
            %   22.5433   -0.3071   -0.0072
            %   -0.3071    1.9685    0.4322
            %   -0.0072    0.4322    0.1916
            
            n = size(theta, 1);
            M = zeros(n);
            for i = 1: n
                dd_theta = zeros(n, 1);
                dd_theta(i) = 1;
                M(:, i) = OpenChain.inverse_dynamics(theta, zeros(n, 1), dd_theta, ...
                                    [0; 0; 0], [0; 0; 0; 0; 0; 0],Mlist, ...
                                    Glist, Slist);
            end
        end
        function JTFtip = force_at_EE(theta, Ftip, Mlist, Glist, Slist)
            % return: joint forces required to create the EE force Ftip
            % EndEffectorForces.m
            % Returns JTFtip: The joint forces and torques required only to create the end-effector force Ftip.
            n = size(theta, 1);
            JTFtip = OpenChain.inverse_dynamics(zeros(n, 1), zeros(n, 1), ...
                                    [0; 0; 0], Ftip, Mlist, Glist, Slist);
        end
        function Tau = force_computed_at_EE(theta, d_theta, e_int, g, ... 
                                        Mlist, Glist, Slist, ...
                                        r_theta, r_d_theta, r_dd_theta, ...
                                        Kp, Ki, Kd)
            % *** CHAPTER 11: ROBOT CONTROL ***
            % ComputedTorque.m
            % Takes thetalist: n-vector of joint variables,
            %       dthetalist: n-vector of joint rates,
            %       eint: n-vector of the time-integral of joint errors,
            %       g: Gravity vector g,
            %       Mlist: List of link frames {i} relative to {i-1} at the home 
            %              position,
            %       Glist: Spatial inertia matrices Gi of the links,
            %       Slist: Screw axes Si of the joints in a space frame, in the format
            %              of a matrix with the screw axes as the columns,
            %       thetalistd: n-vector of reference joint variables,
            %       dthetalistd: n-vector of reference joint velocities,
            %       ddthetalistd: n-vector of reference joint accelerations,
            %       Kp: The feedback proportional gain (identical for each joint),
            %       Ki: The feedback integral gain (identical for each joint),
            %       Kd: The feedback derivative gain (identical for each joint).
            % Returns taulist: The vector of joint forces/torques computed by the 
            %                  feedback linearizing controller at the current instant.
            e = r_theta - theta;
            Tau = OpenChain.mass_matrix(theta, Mlist, Glist, Slist) ...
                * (Kp * e + Ki * (e_int + e) + Kd * (r_d_theta - d_theta)) ...
                + OpenChain.inverse_dynamics(theta, d_theta, r_dd_theta, g, ...
                                            zeros(6,1), Mlist, Glist, Slist);
        end
        function G = gravity_forces(theta, g, Mlist, Glist, Slist)
        % Returns grav: The joint forces/torques required to overcome gravity at 
        %               theta
            n = size(theta, 1);
            G = OpenChain.inverse_dynamics(theta,   zeros(n, 1), ...
                                                    zeros(n, 1), g, ...
                                                    zeros(6, 1), Mlist, Glist, Slist);
        end
        function dd_theta = forward_dynamics(theta, d_theta, tau, g, Ftip, Mlist, Glist, Slist)
            % Returns ddthetalist: The resulting joint accelerations.
            % This function computes ddthetalist by solving:
            % Mlist(thetalist) * ddthetalist = taulist - c(thetalist,dthetalist) ...
            %                                  - g(thetalist) - Jtr(thetalist) * Ftip
            dd_theta = OpenChain.mass_matrix(theta, Mlist, Glist, Slist) \ ( ...
                        tau ...
                        - OpenChain.vel_qualdratic_force(theta, d_theta, Mlist, Glist, Slist) ...
                        - OpenChain.gravity_forces(theta, g, Mlist, Glist, Slist) ...
                        - OpenChain.force_at_EE(theta, Ftip, Mlist, Glist, Slist) ...
                    );
        end
        function traj = trajectory_joint(thetaS, thetaE, Tf, N, method)
            timegap = Tf / (N - 1);
            traj = zeros(size(thetaS, 1), N);
            for i = 1: N
                if method == 3
                    s = CubicTimeScaling(Tf, timegap * (i - 1));
                else
                    s = QuinticTimeScaling(Tf, timegap * (i - 1));
                end
                traj(:, i) = thetaS + s * (thetaE - thetaS);
            end
            traj = traj';
        end
        function taumat = trajectory_ID(theta, d_theta, dd_theta, g, Ftip, Mlist, Glist, Slist)
            thetamat = theta';
            dthetamat = d_theta';
            ddthetamat = dd_theta';
            Ftipmat = Ftip';
            taumat = thetamat;
            for i = 1: size(thetamat, 2)
            taumat(:, i) ...
            = InverseDynamics(thetamat(:, i), dthetamat(:, i), ddthetamat(:, i), ...
                                g, Ftipmat(:, i), Mlist, Glist, Slist);
            end
            taumat = taumat';
        end
    end
end
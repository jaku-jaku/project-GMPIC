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
        function G = body_forward_kinematics(M, Blist, theta)
            % *** CHAPTER 9: TRAJECTORY GENERATION ***
            % Takes M: The home configuration (position and orientation) of the end-
            %          effector,
            %       Blist: The joint screw axes in the end-effector frame when the
            %              manipulator is at the home position,
            %       theta: A list of joint coordinates.
            % Returns G: The end-effector configuration resulting from forward kinematics
            %            (position and orientation) represented in SE(3).

            G = M;
            for i = 1: size(theta, 1)
                G = G * Lie.exp_map_SE3_from_R6(Blist(:, i) * theta(i));
            end
        end
        function G = spatial_forward_kinematics(M, Slist, theta)
            % Takes M: the home configuration (position and orientation) of the 
            %          end-effector,
            %       Slist: The joint screw axes in the space frame when the manipulator
            %              is at the home position,
            %       theta: A list of joint coordinates.
            % Returns T in SE(3) representing the end-effector frame, when the joints 
            % are at the specified coordinates (i.t.o Space Frame).
            G = M;
            for i = size(theta): -1: 1
                G = Lie.exp_map_SE3_from_R6(Slist(:, i) * theta(i)) * G;
            end
        end
        function Js = spatial_jacobian(Slist, theta)
            % *** CHAPTER 5: VELOCITY KINEMATICS AND STATICS ***
            % Takes Slist: The joint screw axes in the space frame when the manipulator
            %              is at the home position, in the format of a matrix with the
            %              screw axes as the columns,
            %       theta: A list of joint coordinates. 
            % Returns the corresponding space Jacobian (6xn real numbers).
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
            % taulist = Mlist(theta) * ddthetalist + c(theta, dthetalist) ...
            %           + g(theta) + Jtr(theta) * Ftip

            n = size(theta, 1);
            Mi = eye(4);
            Ai = zeros(6, n);
            AdTi = zeros(6, 6, n + 1);
            Vi = zeros(6, n + 1);
            Vdi = zeros(6, n + 1);
            Vdi(1:3, 1) = -g;
            AdTi(:, :, n + 1) = Lie.Ad_SE3_from_SE3(Lie.inverse_SE3(Mlist(:, :, n + 1)));
            Fi = Ftip;
            taulist = zeros(n, 1);
            for i=1: n    
                Mi = Mi * Mlist(:, :, i);
                Ai(:, i) = Lie.Ad_SE3_from_SE3(Lie.inverse_SE3(Mi)) * Slist(:, i);    
                % [A generic one]:
                AdTi(:, :, i) = Lie.Ad_SE3_from_SE3( ...
                    Lie.exp_map_SE3_from_R6(Ai(:, i) * -theta(i)) * Lie.inverse_SE3(Mlist(:, :, i)) ...
                );
                % [Assuming, Ai(:,i) contains a unit axis vector: w_R3]:
                % AdTi(:, :, i) = Lie.Ad_SE3_from_SE3( ...
                %     Lie.exp_map_SE3_from_R6xR(Ai(:, i), -theta(i)) * Lie.inverse_SE3(Mlist(:, :, i)) ...
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
            % Returns c: The vector c(theta,dthetalist) of Coriolis and centripetal
            %            terms for a given theta and dthetalist.
            % This function calls InverseDynamics with g = 0, Ftip = 0, and 
            % ddthetalist = 0.
            
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
            % Returns M: The numerical inertia matrix M(theta) of an n-joint serial
            %            chain at the given configuration theta.
            % This function calls InverseDynamics n times, each time passing a 
            % ddthetalist vector with a single element equal to one and all other 
            % inputs set to zero. Each call of InverseDynamics generates a single 
            % column, and these columns are assembled to create the inertia matrix.
            % Example Input (3 Link Robot):
            
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
            % Takes theta: n-vector of joint variables,
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
            % Mlist(theta) * ddthetalist = taulist - c(theta,dthetalist) ...
            %                                  - g(theta) - Jtr(theta) * Ftip
            dd_theta = OpenChain.mass_matrix(theta, Mlist, Glist, Slist) \ ( ...
                        tau ...
                        - OpenChain.vel_qualdratic_force(theta, d_theta, Mlist, Glist, Slist) ...
                        - OpenChain.gravity_forces(theta, g, Mlist, Glist, Slist) ...
                        - OpenChain.force_at_EE(theta, Ftip, Mlist, Glist, Slist) ...
                    );
        end
        function [theta, success] = spatial_inverse_kinematics(Slist, M, T, theta_0, e_w, e_v)
            % *** CHAPTER 6: INVERSE KINEMATICS ***
            % Takes Slist: The joint screw axes in the space frame when the manipulator
            %              is at the home position, in the format of a matrix with the
            %              screw axes as the columns,
            %       M: The home configuration of the end-effector,
            %       T: The desired end-effector configuration Tsd,
            %       theta_0: An initial guess of joint angles that are close to 
            %                   satisfying Tsd,
            %       e_w: A small positive tolerance on the end-effector orientation 
            %             error. The returned joint angles must give an end-effector 
            %             orientation error less than e_w,
            %       e_v: A small positive tolerance on the end-effector linear position 
            %           error. The returned joint angles must give an end-effector 
            %           position error less than e_v.
            % Returns theta: Joint angles that achieve T within the specified 
            %                    tolerances,
            %         success: A logical value where TRUE means that the function found
            %                  a solution and FALSE means that it ran through the set 
            %                  number of maximum iterations without finding a solution
            %                  within the tolerances e_w and e_v.
            % Uses an iterative Newton-Raphson root-finding method.
            % The maximum number of iterations before the algorithm is terminated has 
            % been hardcoded in as a variable called maxiterations. It is set to 20 at 
            % the start of the function, but can be changed if needed.  
            theta = theta_0;
            i = 0;
            maxiterations = 20;
            Tsb = OpenChain.spatial_forward_kinematics(M, Slist, theta);
            Tsb_se3 = Lie.log_se3_from_SE3(Lie.inverse_SE3(Tsb) * T);
            Vs = Lie.Ad_SE3_from_SE3(Tsb) * Lie.vee_R6_from_se3(Tsb_se3)';
            err = norm(Vs(1: 3)) > e_v || norm(Vs(4: 6)) > e_w;
            while err && i < maxiterations
                theta = theta + pinv(OpenChain.spatial_jacobian(Slist, theta)) * Vs;
                i = i + 1;
                Tsb = OpenChain.spatial_forward_kinematics(M, Slist, theta);
                Tsb_se3 = Lie.log_se3_from_SE3(Lie.inverse_SE3(Tsb) * T);
                Vs = Lie.Ad_SE3_from_SE3(Tsb) * Lie.vee_R6_from_se3(Tsb_se3)';
                err = norm(Vs(1: 3)) > e_v || norm(Vs(4: 6)) > e_w;
            end
            success = ~ err;
        end
        function [theta, success] = body_inverse_kinematics(Blist, M, T, theta_0, e_w, e_v)
            % Takes Blist: The joint screw axes in the end-effector frame when the
            %              manipulator is at the home position, in the format of a 
            %              matrix with the screw axes as the columns,
            %       M: The home configuration of the end-effector,
            %       T: The desired end-effector configuration Tsd,
            %       thetalist0: An initial guess of joint angles that are close to 
            %                   satisfying Tsd,
            %       e_w: A small positive tolerance on the end-effector orientation
            %             error. The returned joint angles must give an end-effector 
            %             orientation error less than eomg,
            %       e_v: A small positive tolerance on the end-effector linear position 
            %           error. The returned joint angles must give an end-effector
            %           position error less than ev.
            % Returns theta: Joint angles that achieve T within the specified 
            %                    tolerances,
            %         success: A logical value where TRUE means that the function found
            %                  a solution and FALSE means that it ran through the set 
            %                  number of maximum iterations without finding a solution
            %                  within the tolerances eomg and ev.
            % Uses an iterative Newton-Raphson root-finding method.
            % The maximum number of iterations before the algorithm is terminated has 
            % been hardcoded in as a variable called maxiterations. It is set to 20 at 
            % the start of the function, but can be changed if needed.  
            
            theta = theta_0;
            i = 0;
            maxiterations = 20;
            Tb = OpenChain.body_forward_kinematics(M, Blist, theta);
            Vb = Lie.vee_R6_from_se3(Lie.log_se3_from_SE3(Lie.inverse_SE3(Tb) * T))';
            
            err = norm(Vb(1: 3)) > e_v || norm(Vb(4: 6)) > e_w;
            while err && i < maxiterations
                theta = theta + pinv(OpenChain.body_jacobian(Blist, theta)) * Vb;
                i = i + 1;
                Tb = OpenChain.body_forward_kinematics(M, Blist, theta);
                Vb = Lie.vee_R6_from_se3(Lie.log_se3_from_SE3(Lie.inverse_SE3(Tb) * T))';
                err = norm(Vb(1: 3)) > e_v || norm(Vb(4: 6)) > e_w;
            end
            success = ~ err;
        end
        function s = time_scaling(Tf, t, mode)
            switch mode
                case "linear"
                    s = t / Tf;
                case "cubic" 
                    % 3rd poly motion that begins and ends at zero velocity.
                    s = 3 * (t / Tf) ^ 2 - 2 * (t / Tf) ^ 3;
                case "quintic" 
                    % 5th poly motion that begins and ends at zero velocity and zero acceleration.
                    s = 10 * (t / Tf) ^ 3 - 15 * (t / Tf) ^ 4 + 6 * (t / Tf) ^ 5;
            end
        end
        function [theta_1, d_theta_1] = euler_stepping(theta, d_theta, dd_theta, dt)
            theta_1 = theta + d_theta * dt;
            d_theta_1 = d_theta + dd_theta * dt;
        end
    end
end
classdef WAM_model
    %% USAGE: pre-req:  c0 = WAM_model.model_init_zero_config(); 
    methods(Static)
        %% Dynamic Sim Function:
        function data = dynamic_sim(c0, INIT_ANGLES, LIST_OF_JOINT_ANGLES, LIST_OF_F_EE, dT, FB_Ctrl_CONFIG)
            % Example Params:
            %       - FB_Ctrl_CONFIG = struct("Kp", 1, "Ki", 0.8, "Kd", 1);
            [N_angles,N_sim] = size(LIST_OF_JOINT_ANGLES);
            WAM_model.assert_if_jnts(c0,N_angles);
            % set initial angles:
            % c0 = c0.set_angle(INIT_ANGLES);
            data = struct( ...
                "theta",zeros(N_angles,N_sim+1), ...
                "d_theta",zeros(N_angles,N_sim+1), ...
                "dd_theta",zeros(N_angles,N_sim+1) ...
            );
            data.theta(:,1) = INIT_ANGLES;
            d_theta_r_k_prev = 0;
            dd_theta_r_k_prev = 0;

            helper.textprogressbar('Stepping through simulations: ');
            % stepping through:
            for k=1:N_sim
                % --- 
                % 1. Stepping:
                % ---
                % 1.0 current state:
                theta_k = data.theta(:,k);
                d_theta_k = data.d_theta(:,k);

                % 1.1 joint reference:
                theta_r_k = LIST_OF_JOINT_ANGLES(:,k);
                theta_r_k = WAM_model.filter_angles(c0, theta_r_k);
                theta_r_k = WAM_model.joint_constraints(c0, theta_r_k);

                % 1.2 compute delta reference:
                d_theta_r_k = (theta_r_k - d_theta_r_k_prev)/dT;
                dd_theta_r_k = (d_theta_r_k - dd_theta_r_k_prev)/dT;
                
                % --- 
                % 2. Controller:
                % ---
                % 2.1 TODO: these disturbance can be injected somehow
                % force_at_jnt = zeros(N_angles, 1); % [NOTE: zero external force]d
                % force_at_EE = zeros(6, 1);  % [NOTE: zero external force] : [Fx; Fy; Fz; tx; ty; tz];
                force_at_EE = LIST_OF_F_EE(:,k);

                % 2.2 feedback linearization controller:
                e_int = (theta_r_k - theta_k) * dT; % linear approx (1-Step Euler)
                T_jnt = WAM_model.feedback_linearization_controller( ...
                    c0, theta_k, d_theta_k, e_int,...
                    theta_r_k, d_theta_r_k, dd_theta_r_k,...
                    FB_Ctrl_CONFIG.Kp, FB_Ctrl_CONFIG.Ki, FB_Ctrl_CONFIG.Kd);
                % --- 
                % 2.b MPC Controller:
                % ---
                % compute Task Frame:
                Xd_SE3 = WAM_model.compute_spatial_FK(c0, theta_r_k);
                X00_SE3 = WAM_model.compute_spatial_FK(c0, theta_k);
                
                % compute error state MPC
                % u = WAM_model.eMPC_casadi(xi0, xid, X00_SE3, Xd_SE3, param)
                
                % --- 
                % 3. Forward sim:
                % ---
                % 3.2 compute dynamics: 
                dd_theta_k = WAM_model.forward_dynamics(c0, theta_k, d_theta_k, force_at_EE, T_jnt);

                % 3.3 integrate to next time step (Approx.): (1-Step Euler)
                % TODO: Alternative Integrations: n-Euler, Ode15, RK 
                d_theta_k_next = d_theta_k + dd_theta_k * dT;
                theta_k_next = theta_k + d_theta_k_next * dT;
                
                % joint constraints:
                theta_k_next = WAM_model.filter_angles(c0, theta_k_next);
                theta_k_next = WAM_model.joint_constraints(c0, theta_k_next);

                % update:
                data.theta(:,k+1) = theta_k_next;
                data.d_theta(:,k+1) = d_theta_k_next;
                data.dd_theta(:,k) = dd_theta_k;

                % cache for next itr:
                d_theta_r_k_prev = theta_r_k;
                dd_theta_r_k_prev = d_theta_r_k;
                
                % kick the dog:
                helper.textprogressbar(k/N_sim);
            end
            helper.textprogressbar('done');
        end
        %% MPC:
        function u = eMPC_casadi(xi0, xid, X00, Xd, param)
            % Xd = [quat2rotm(xd(1:4)), xd(5:7)';...
            %     0,0,0,1];
            % X00 = [quat2rotm(x0(1:4)'), x0(5:7);...
            %     0,0,0,1];
            
            % let xi0 = x0(8:end) \in se3 coord, and X00, Xd \in SE3:
            
            % The inverse of Left error defined in paper is equivalent.
            X00 =  X00^(-1) * Xd;
            X0 = logm(X00);
            
            p0 = [X0(3,2); X0(1,3); X0(2,1); X0(1:3,4);xi0(8:end)]; %% initial error
            %%
            I = param.I;
            dt = param.dt;
            Nx = param.Nx;
            %% use Casadi:
            helper.includeCasadi();
            import casadi.*
            %%
            x = SX.sym('x', [param.Nx, param.Nt + 1]);
            u = SX.sym('u', [param.Nu, param.Nt]);
            
            g_dynamics = [];
            g_dynamics_min = [];
            g_dynamics_max = [];
            J = 0;
            for k = 1:param.Nt
                % xi_bar = I * xid(k, :)';
                xi_bar = I * xi0;
                G = ad_se3_from_R6(xi_bar); % (v,w)
                
                coadjoint = Lie.ad_se3_from_R6(xi0)';
                H = I^-1 *(coadjoint * I + G);
                
                Ac = [-Lie.ad_se3_from_R6(xid(k, :)),   -eye(6);...
                    zeros(6), H];
                Bc = [zeros(6);...
                    I^-1];
                
                b = -I^(-1) * G * xi0;
                hc = [xid(k,:)'; b];
                
                Ad = expm(Ac * dt);
                Bd = (eye(param.Nx) * dt + Ac * dt^2 / 2 + Ac^2 * dt^3 / 6) * Bc;
                hd = (eye(param.Nx) * dt + Ac * dt^2 / 2 + Ac^2 * dt^3 / 6) * hc;
                
                X_k = x(:, k);
                X_kp1 = x(:, k+1);
                u_k = u(:, k);
                
                g_dynamics = [g_dynamics;...
                    X_kp1 - Ad * X_k - Bd * u_k - hd];
                g_dynamics_min = [g_dynamics_min;...
                    zeros(Nx,1)];
                g_dynamics_max = [g_dynamics_max;...
                    zeros(Nx,1)];
                
                C = eye(12);
                C(7:12, 1:6) = -Lie.ad_se3_from_R6(xid(k,:));
                d = zeros(12,1);
                d(7:12) = xid(k,:)';
                J = J + 0.5 * (C * x(:,k) - d)' * param.Q * (C * x(:,k) - d) + 0.5 * u(:,k)' * param.R * u(:,k);
            end
            
            P = idare(Ad, Bd, param.Q,param.R,[], []);
            J = J + 0.5 * (C * X_kp1 - d)' * P * (C * X_kp1 - d);
            
            g_boundary = [];
            g_boundary = [g_boundary; x(:,1) - p0]; % initial error
            
            g = [g_dynamics; g_boundary];
            g_min = [g_dynamics_min; zeros(param.Nx, 1)];
            g_max = [g_dynamics_max; zeros(param.Nx, 1)];
            
            b_input_min = repmat(param.umin, [param.Nt, 1]);
            b_input_max = repmat(param.umax, [param.Nt, 1]);
            
            x_max = [b_input_max; inf(param.Nx * (param.Nt+1), 1)];
            x_min = [b_input_min; -inf(param.Nx * (param.Nt+1), 1)];
            
            input = [u(:); x(:)];
            %%
            nlp = struct('x',input, 'f',J, 'g', g);
            % help nlpsol
            S = nlpsol('S', 'ipopt', nlp);
            
            x0 = randn(size(input));
            r = S('x0',x0,'lbg',g_min,'ubg',g_max,'lbx',x_min,'ubx',x_max);
            sol = full(r.x);
            u = sol(1:6);
        end
        %% Configuration dependent:
        function T_jnt = feedback_linearization_controller(c0, ...
                theta, d_theta, e_int, r_theta, r_d_theta, r_dd_theta, Kp, Ki, Kd)
            e = r_theta - theta;
            e_d = r_d_theta - d_theta;
            m_mass_MR = OpenChainMR.mass_matrix(theta, c0.Mlist, c0.Glist, c0.Slist);
            % feedback linearing controller
            T_jnt_fb = OpenChainMR.inverse_dynamics(theta, d_theta, r_dd_theta, ...
                              c0.g, zeros(6,1), c0.Mlist, c0.Glist, c0.Slist);
            % inverse dynamics controller:
            %       = feedforward + feedback linearing controller
            T_jnt = m_mass_MR * (Kp * e + Ki * (e_int + e) + Kd * e_d) + T_jnt_fb;
            % equivalent to:
            % T_jnt = OpenChainMR.computed_torque_at_joints(theta, d_theta, e_int, ... 
            %                             c0.g, c0.Mlist, c0.Glist, c0.Slist, ...
            %                             r_theta, r_d_theta, r_dd_theta, ...
            %                             Kp, Ki, Kd);
        end
        function dd_theta = forward_dynamics(c0, theta, d_theta, F_EE_twist, T_jnt)
            % compute config-dependent dyn variables:
            m_mass_MR = OpenChainMR.mass_matrix(theta, c0.Mlist, c0.Glist, c0.Slist);
            m_gravity_MR = OpenChainMR.gravity_forces(theta, c0.g, c0.Mlist, c0.Glist, c0.Slist);
            m_Coriolis_MR = OpenChainMR.vel_qualdratic_force(theta, d_theta, c0.Mlist, c0.Glist, c0.Slist);
            m_EE_force_MR = OpenChainMR.EE_Force_to_Joint_Torques(theta, F_EE_twist, c0.Mlist, c0.Glist, c0.Slist);
            % compute dynamics:
            dd_theta = m_mass_MR \ (T_jnt - m_Coriolis_MR - m_gravity_MR - m_EE_force_MR);
        end
        %% Stationary Functions:
        function plot_zero_config(c0)
            WAM_model.plot_at(c0, zeros(c0.N_JNTS,1))
        end
        function theta = filter_angles(c0, JOINT_ANGLES)
            % apply safe filter:
            theta = c0.VALID_JNTS_FILTER_ * JOINT_ANGLES;
        end
        function theta = joint_constraints(c0, theta)
            % TODO: apply joint constraints:
            theta = max(c0.JOINT_LIMITS_(:,1),theta);
            theta = min(c0.JOINT_LIMITS_(:,2),theta);
        end
        function plot_at(c0, JOINT_ANGLES)
            % --- 
            % FK SIM:
            % --- 
            G_SE3_wam_spatial_ours_ = WAM_model.compute_spatial_FK_for_all_joints(c0, JOINT_ANGLES);
            WAM_model.plot_WAM(true, G_SE3_wam_spatial_ours_, c0.G_SE3_s_0);
        end
        function assert_if_jnts(c0,N_jnts)
            assert(N_jnts==c0.N_JNTS, ...
                sprintf("Joint Angles Array should be in form of %dxN_STEPS",c0.N_JNTS))
        end
        function G_SE3_EE = compute_spatial_FK(c0, theta)
            % - Ours:
            % exp_xi_theta_in_SE3_ = OpenChain.batch_screw_SE3_from_twist_angle_R6xR(c0.xi_R6_s_, theta);
            % G_SE3_EE_spatial_traj_{k} = OpenChain.compute_Spatial_FK(exp_xi_theta_in_SE3_, c0.G_SE3_s_{end});
            % - MR:
            G_SE3_EE = OpenChainMR.spatial_forward_kinematics(c0.G_SE3_s_{end}, c0.Slist, theta);
        end
        function G_SE3_EE_ = compute_spatial_FK_for_all_joints(c0, theta)
            % - Ours:
            exp_xi_theta_in_SE3_ = OpenChain.batch_screw_SE3_from_twist_angle_R6xR(c0.xi_R6_s_, theta);
            G_SE3_EE_ = OpenChain.compute_Spatial_FK_for_all_joints(exp_xi_theta_in_SE3_, c0.G_SE3_s_);
        end
        function G_SE3_EE_spatial_traj_ = compute_EE_trajectory_with(c0, list_of_thetas)
            [N_jnts, N_rec] = size(list_of_thetas);
            WAM_model.assert_if_jnts(c0,N_jnts);
            G_SE3_EE_spatial_traj_ = cell(N_jnts,1);
            for k=1:N_rec
                theta = list_of_thetas(:,k);
                G_SE3_EE_spatial_traj_{k} = WAM_model.compute_spatial_FK(c0, theta);
            end
        end
        function animate_with(c0,list_of_thetas,VIEW_DIMENSION,VIEW_SIZE,VIEW_ANGLE,DIR,TAG,FPS)
            filename = sprintf("WAM_model_animated_%s",TAG);
            helper.createRecorder(DIR,filename,100,FPS);
            [N_jnts, N_rec] = size(list_of_thetas);
            WAM_model.assert_if_jnts(c0,N_jnts);
            helper.newFigure(-1);
            title(filename);
            helper.textprogressbar('generating videos: ');

            % - plot and record:
            for i = 1:N_rec
                clf("reset")
                helper.setPlotSize(VIEW_DIMENSION);
                view(VIEW_ANGLE)
                utils.plot3DBoundary(VIEW_SIZE)
                axis equal;
                hold on;
                WAM_model.plot_at(c0, list_of_thetas(:,i));
                grid on;
                hold off;
                drawnow;
                helper.recordRecorder();
                helper.textprogressbar(i/N_rec);
                pause(0.000001);
            end
            helper.textprogressbar('done');
            helper.loginfo(sprintf("Done Recording @ %s",filename))
            helper.terminateRecorder();
        end
        %% [Config independent] Model Initialization:
        function c0 = model_init_zero_config()
            c0 = struct();
            % --- 
            % Initialization:
            % --- 
            [   c0.M_R6x6_, c0.G_SE3_, ...
                c0.xi_R6_s_, c0.G_SE3_s_, c0.xi_R6_s_0, c0.G_SE3_s_0, ... % spatial
                c0.xi_R6_b_, c0.G_SE3_b_, c0.xi_R6_b_0, c0.G_SE3_b_0, ... % body
                c0.JOINT_LIMITS_, c0.VALID_JNTS_FILTER_ ...
            ] = WAM_model.model_initialization();
            c0.N_JNTS = length(c0.xi_R6_s_);
            
            % convert to MR:
            c0.Slist = cat(2, c0.xi_R6_s_{:});
            c0.Mlist = cat(3, c0.G_SE3_s_0, c0.G_SE3_{:});
            c0.Glist = cat(3, c0.M_R6x6_{:});
            c0.Blist = cat(2, c0.xi_R6_b_{:});

            % other conatnt
            c0.g = common.GRAVITY; 
        end
        function [M_R6x6_, G_SE3_, ...
                xi_R6_s_, G_SE3_s_, xi_R6_s_0, G_SE3_s_0, ... % spatial
                xi_R6_b_, G_SE3_b_, xi_R6_b_0, G_SE3_b_0, ... % body
                joint_limits_, valid_jnts_filter] = model_initialization()
            %% [Forward for Spatial]:
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
            joint_limits_ = zeros(N_links,2); 
            G_SE3_ = cell(1,N_links); 
            % world_coord-->0 (wam base):
            w_joint_axis_0 = [0;0;1]; 
            G_SE3_s_0 = WAM_Spatial_0;
            w_R3_s_0 = G_SE3_s_0(1:3,1:3) * w_joint_axis_0;
            xi_R6_s_0 = RBT.twist_R6_from_axis_point(w_R3_s_0, G_SE3_s_0(1:3,4));
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
                joint_limits_(i,:) = common.WAM(i).tip_frame_revolute_joint.limits;
            
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
                    G_SE3_s_{i} = G_SE3_s_0 * G_SE3_{i}; % <--- initial wam base joint tip frame (J1)
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
            
            %% [iterating links backward]:
            G_SE3_b_ = cell(1,N_links);
            xi_R6_b_ = cell(1,N_links);
            
            G_TOOL_0 = eye(4); % the end frame is the body/tool frame (the identity)
            % todo: please make sure the produced link can provide same jacobian
            for i = N_links:-1:1
                helper.loginfo(sprintf("> Indexing Links @ %d-->%d",i,i-1));
                % - propagate:
                if i == N_links
                    G_SE3_b_{i} = G_TOOL_0;
                else
                    G_SE3_b_{i} = G_SE3_b_{i+1} * RBT.inverse_SE3(G_SE3_{i+1}); 
                end
            
                % - obtain global representation of the joint rotation axis
                w_joint_axis_i = common.WAM(i).tip_frame_revolute_joint.axis';
                
                % - twist:
                w_R3_b_i = G_SE3_b_{i}(1:3,1:3) * w_joint_axis_i; 
                xi_R6_b_{i} = RBT.twist_R6_from_axis_point(w_R3_b_i, G_SE3_b_{i}(1:3,4));
            end
            G_SE3_b_0 = G_SE3_b_{1} * RBT.inverse_SE3(G_SE3_{1});
            w_R3_b_0 = G_SE3_b_0(1:3,1:3) * w_joint_axis_0; 
            xi_R6_b_0 = RBT.twist_R6_from_axis_point(w_R3_b_0, G_SE3_b_0(1:3,4));
            xi_R6_b = cat(2,xi_R6_b_{:})
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
            
%             grid on;
%             axis equal;
            xlabel('X Axis')
            ylabel('Y Axis')
            zlabel('Z Axis')
        end
    end
end

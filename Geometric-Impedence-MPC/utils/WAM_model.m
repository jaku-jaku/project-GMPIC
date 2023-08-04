classdef WAM_model
    properties(Constant)
        GRAVITY = [0, 0, 9.81]';
    end
    properties
        % init values:
        WAM_Spatial_0
        xi_R6_s_
        G_SE3_s_
        G_SE3_
        M_R6x6_
        JOINT_LIMITS_
        VALID_JNTS_FILTER_
        N_JNTS_
        % [MR]
        Slist 
        Mlist
        Glist
        % dynamic values:
%         m_J_wam_spatial_MR
%         m_mass_MR
%         m_gravity_MR
%         m_Coriolis_MR
        m_theta
%         m_d_theta
%         m_dd_theta
    end
    methods
        function obj = WAM_model()
            % --- 
            % Initialization:
            % --- 
            [obj.WAM_Spatial_0, obj.xi_R6_s_, obj.G_SE3_s_, obj.M_R6x6_, ...
                obj.G_SE3_, obj.JOINT_LIMITS_, obj.VALID_JNTS_FILTER_ ...
            ] = WAM_model.model_initialization();
            obj.N_JNTS_ = length(obj.xi_R6_s_);
            
            % convert to MR:
            obj.Slist = cat(2, obj.xi_R6_s_{:});
            obj.Mlist = cat(3, eye(4), obj.G_SE3_{:});
            obj.Glist = cat(3, obj.M_R6x6_{:});

            % init:
            obj.m_theta = zeros(obj.N_JNTS_,1);
        end
        %% Time Dependent Functions:
        function [obj,data] = dynamic_sim(obj, INIT_ANGLES, LIST_OF_JOINT_ANGLES)
            [N_angles,N_sim] = size(LIST_OF_JOINT_ANGLES);
            obj.assert_if_jnts(N_jnts);
            % set initial angles:
            obj = obj.set_angles(INIT_ANGLES);
            data = struct( ...
                "theta",zeros(N_angles,N_sim+1), ...
                "d_theta",zeros(N_angles,N_sim+1), ...
                "dd_theta",zeros(N_angles,N_sim+1) ...
            );
            data.theta(:,1) = INIT_ANGLES;
            % stepping through:
            for k=1:N_sim
                % 1. joint request:
                theta_r_k = LIST_OF_JOINT_ANGLES(:,k);
                theta_r_k = obj.filter_angles(theta_r_k);
                theta_r_k = obj.joint_constraints(theta_r_k);
                % 2. compute delta:
                d_theta_r_k = theta_r_k - obj.m_theta;
                
                % 3. compute configuration dependent var:
                % --- 
                % Jacobi:
                % --- 
                obj.m_J_wam_spatial_MR = OpenChainMR.spatial_jacobian(obj.Slist, theta_r_k);
                
                % --- 
                % dynamics:
                % --- 
                m_mass_MR = OpenChainMR.mass_matrix(theta_r_k, obj.Mlist, obj.Glist, obj.Slist);
                m_gravity_MR = OpenChainMR.gravity_forces(theta_r_k, obj.GRAVITY, obj.Mlist, obj.Glist, obj.Slist);
                m_Coriolis_MR = OpenChainMR.vel_qualdratic_force(theta_r_k, d_theta_r_k, obj.Mlist, obj.Glist, obj.Slist);

                % compute dynamics:
                % TODO: these disturbance can be injected somehow
                force_at_jnt = zeros(N_angles, 1); % [NOTE: zero external force]
                force_at_EE = zeros(N_angles, 1);  % [NOTE: zero external force]
                dd_theta_k = m_mass_MR \ (force_at_jnt - m_Coriolis_MR - m_gravity_MR - force_at_EE);

                % update:
                data.theta(:,1) = theta_r_k;
                data.d_theta(:,1) = d_theta_r_k;
                data.dd_theta(:,1) = dd_theta_k;
            end
        end
        %% Set Functions:
        function obj = set_angle(obj, JOINT_ANGLES)
            obj.m_theta = obj.filter_angles(JOINT_ANGLES);
            obj.m_theta = obj.joint_constraints(obj.m_theta);
        end
        %% Stationary Functions:
        function theta = filter_angles(obj, JOINT_ANGLES)
            % apply safe filter:
            theta = obj.VALID_JNTS_FILTER_ * JOINT_ANGLES;
        end
        function theta = joint_constraints(obj, theta)
            % TODO: apply joint constraints:
            theta = max(obj.JOINT_LIMITS_(:,1),theta);
            theta = min(obj.JOINT_LIMITS_(:,2),theta);
        end
        function plot_now(obj)
            obj.plot_at(obj.m_theta)
        end
        function plot_at(obj, JOINT_ANGLES)
            % --- 
            % FK SIM:
            % --- 
            exp_xi_theta_in_SE3_ = OpenChain.batch_screw_SE3_from_twist_angle_R6xR(obj.xi_R6_s_, JOINT_ANGLES);
            G_SE3_wam_spatial_ours_ = OpenChain.compute_Spatial_FK_for_all_joints(exp_xi_theta_in_SE3_, obj.G_SE3_s_);
            WAM_model.plot_WAM(true, G_SE3_wam_spatial_ours_, obj.WAM_Spatial_0);
        end
        function assert_if_jnts(obj,N_jnts)
            assert(N_jnts==obj.N_JNTS_, ...
                sprintf("Joint Angles Array should be in form of %dxN_STEPS",obj.N_JNTS_))
        end
        function G_SE3_EE_spatial_traj_ = compute_EE_trajectory_with(obj, list_of_thetas)
            [N_jnts, N_rec] = size(list_of_thetas);
            obj.assert_if_jnts(N_jnts);
            G_SE3_EE_spatial_traj_ = cell(N_jnts,1);
            for k=1:N_rec
                theta = list_of_thetas(:,k);
                exp_xi_theta_in_SE3_ = OpenChain.batch_screw_SE3_from_twist_angle_R6xR(obj.xi_R6_s_, theta);
                G_SE3_EE_spatial_traj_{k} = OpenChain.compute_Spatial_FK(exp_xi_theta_in_SE3_, obj.G_SE3_s_{end});
            end
        end
        function animate_with(obj,list_of_thetas,VIEW_SIZE,VIEW_ANGLE,DIR,TAG,FPS)
            filename = sprintf("WAM_model_%s",TAG);
            helper.createRecorder(DIR,filename,100,FPS);
            [N_jnts, N_rec] = size(list_of_thetas);
            obj.assert_if_jnts(N_jnts);
            helper.newFigure(-1);
            title(filename);
            helper.textprogressbar('generating videos: ');
            
            % plot and record:
            for i = 1:N_rec
                clf("reset")
                view(VIEW_ANGLE)
                utils.plot3DBoundary(VIEW_SIZE)
                axis equal;
                hold on;
                obj.plot_at(list_of_thetas(:,i));
                grid on;
                hold off;
                drawnow;
                helper.recordRecorder();
                helper.textprogressbar(i);
                pause(0.000001);
            end
            helper.textprogressbar('done');
            helper.loginfo(sprintf("Done Recording @ %s",filename))
            helper.terminateRecorder();
            
            helper.textprogressbar('terminated');
        end
    end
    %% Here are the model independent functions:
    methods(Static)
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
            joint_limits_ = zeros(N_links,2); 
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
            
%             grid on;
%             axis equal;
            xlabel('X Axis')
            ylabel('Y Axis')
            zlabel('Z Axis')
        end
    end
end

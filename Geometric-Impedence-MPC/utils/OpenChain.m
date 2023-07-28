classdef OpenChain
    methods(Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Open Chain Formulation Based on Murray : Mathematical Intro to Robotic Manipulation
        %% --> a more efficient formulation, assuming w and xi is a unit vector~!
        function exp_xi_theta_in_SE3_1_ = batch_screw_SE3_from_twist_angle_R6xR(xi_R6_1_, theta_)
            N_jnts = length(theta_);
            exp_xi_theta_in_SE3_1_ = cell(1,N_jnts);
            for i=1:N_jnts
                % -> compute transformation from twist coordinates and angles applied:
                % |---- (more general):
%                 exp_xi_theta_in_SE3_1_{i} = Lie.exp_map_SE3_from_R6xR(xi_R6_1_{i}, theta_(i));
                % |---- (more efficient):
                exp_xi_theta_in_SE3_1_{i} = RBT.screw_SE3_from_xi_theta(xi_R6_1_{i}, theta_(i));
%                 exp_xi_theta_in_SE3_1_{i} = RBT.screw_SE3_from_axis_point_theta(w_R3_0_{i}, q_R3_0_{i}, theta_(i));
            end
        end
        function G_st_SE3_1_ = compute_FWD_Exponential_Chain(exp_xi_theta_in_SE3_)
            N_jnts = length(exp_xi_theta_in_SE3_);
            G_st_SE3_1_ = cell(1, N_jnts); % -> (placeholder)
            G_st_SE3_1_{1} = exp_xi_theta_in_SE3_{1};
            for i=2:N_jnts
                G_st_SE3_1_{i} = G_st_SE3_1_{i-1} * exp_xi_theta_in_SE3_{i};
            end
        end
        function G_st_SE3_n_ = compute_RVR_Exponential_Chain(exp_xi_theta_in_SE3_)
            N_jnts = length(exp_xi_theta_in_SE3_);
            G_st_SE3_n_ = cell(1, N_jnts); % -> (placeholder)
            G_st_SE3_n_{N_jnts} = exp_xi_theta_in_SE3_{N_jnts};
            for i=N_jnts-1:-1:1
                G_st_SE3_n_{i} = exp_xi_theta_in_SE3_{i} * G_st_SE3_n_{i+1};
            end
        end
        %% % [ Cell Based Batch Process ] 
        function cell_ = cell_prod(A,B)
            [ AN, AM ] = size(A);
            cell_ = cell(AN,AM);
            for i = 1:AN
                for j = 1:AM
                    cell_{i,j} = A{i,j} * B{i,j};
                end
            end
        end
        %% % [ Spatial ] %%% %%% %%% %%% %%% %%% %%%
        function G_s_st_SE3 = compute_Spatial_FK(exp_xi_theta_in_SE3_1_, g0_st_SE3)
            % @G0_SE3_st_EE: zero/home configuration of the ToolFrame: g_st(0)
            N_jnts = length(exp_xi_theta_in_SE3_1_);
            % > g_st(t) = exp_xi_1(t_1) * exp_xi_2(t_2) * ... * exp_xi_n(t_n) * g_st(0)
            G_s_st_SE3 = eye(4); 
            for i=1:N_jnts
                G_s_st_SE3 = G_s_st_SE3 * exp_xi_theta_in_SE3_1_{i};
            end
            G_s_st_SE3 = G_s_st_SE3 * g0_st_SE3; % apply zero config at EE joint
        end
        function G_s_st_SE3_ = compute_Spatial_FK_for_all_joints(exp_xi_theta_in_SE3_1_, g0_st_SE3_)
%             % [ NEW IMPLICIT IMPLEMENTATION ]
%             % @g0_st_SE3_: zero/home configuration for a series of joints: g_st_i(0)
%             % > g_st(t) = exp_xi_1(t_1) * exp_xi_2(t_2) * ... * exp_xi_n(t_n) * g_st(0)
%             g_st_SE3 = OpenChain.compute_FWD_Exponential_Chain(exp_xi_theta_in_SE3_1_);
%             G_s_st_SE3_ = OpenChain.cell_prod(g_st_SE3, g0_st_SE3_);
            
            % [ ORIGINAL EXPLICIT IMPLEMENTATION ]
            % @g0_st_SE3_: zero/home configuration for a series of joints: g_st_i(0)
            N_jnts = length(exp_xi_theta_in_SE3_1_);
            validate.if_dimension("g0_st_SE3_", g0_st_SE3_, [1 N_jnts]);
            G_s_st_SE3_ = cell(1, N_jnts); % -> (placeholder)
            % > g_st(t) = exp_xi_1(t_1) * exp_xi_2(t_2) * ... * exp_xi_n(t_n) * g_st(0)
            g_st_SE3 = eye(4); 
            
            for i=1:N_jnts
                g_st_SE3 = g_st_SE3 * exp_xi_theta_in_SE3_1_{i};
                G_s_st_SE3_{i} = g_st_SE3 * g0_st_SE3_{i}; % apply zero config at each joint
            end
        end
        function J_s_st = compute_Spatial_Jacobian(xi_R6_1_, exp_xi_theta_in_SE3_1_)
            N_jnts = length(xi_R6_1_);
            
            % J = [xi_1, xi_2', ... , xi_N' ]:
            J_s_st = zeros(6,N_jnts);
            
            % [xi_1]:
            J_s_st(:,1) = xi_R6_1_{1};
            g_st_SE3_1_i = eye(4); % init
            
            % [xi_2' ... xi_N']:
            for i=2:N_jnts
                % -> compute current jacobian:
                % compute transformation for 1 ... i-1:
                g_st_SE3_1_i = g_st_SE3_1_i * exp_xi_theta_in_SE3_1_{i-1}; 
                % compute adjoint for 1 ... i-1:
                Ad_1_i = Lie.Ad_SE3_from_SE3(g_st_SE3_1_i);
                % compute jacobian @ xi_i:
                J_s_st(:,i) =  Ad_1_i * xi_R6_1_{i};
            end
        end
        %% % [ Body ] %%% %%% %%% %%% %%% %%% %%%
        function G_b_st_SE3 = compute_Body_FK(exp_xi_theta_in_SE3_body_, g0_st_SE3_body_)
            % @G0_SE3_st_EE: zero/home configuration of the ToolFrame: g_st(0)
            N_jnts = length(exp_xi_theta_in_SE3_body_);
            % > g_st(t) = exp_xi_1(t_1) * exp_xi_2(t_2) * ... * exp_xi_n(t_n) * g_st(0)
            G_b_st_SE3 = g0_st_SE3_body_; 
            for i=1:N_jnts
                G_b_st_SE3 = G_b_st_SE3 * exp_xi_theta_in_SE3_body_{i};
            end
        end
        function G_s_st_SE3_ = compute_Body_FK_for_all_joints(exp_xi_theta_in_SE3_body_, g0_st_SE3_body_)
            % [ NEW IMPLICIT IMPLEMENTATION ]
            % @g0_st_SE3_: zero/home configuration for a series of joints: g_st_i(0)
            % > g_st(t) = exp_xi_1(t_1) * exp_xi_2(t_2) * ... * exp_xi_n(t_n) * g_st(0)
            g_st_SE3 = OpenChain.compute_RVR_Exponential_Chain(exp_xi_theta_in_SE3_body_);
            G_s_st_SE3_ = OpenChain.cell_prod(g0_st_SE3_body_, g_st_SE3);
        end   
        function J_b_st = compute_Body_Jacobian(xi_R6_body_, neg_exp_xi_theta_in_SE3_body_)
            %{
                Expecting:
                    neg_exp_xi_theta_in_SE3_body_ = OpenChain.batch_screw_SE3_from_twist_angle_R6xR(xi_R6_b_, - JOINT_ANGLE);
            %}
            N_jnts = length(xi_R6_body_);
            % J = [xi_1#, xi_2#, ... , xi_N# ]:
            J_b_st = zeros(6,N_jnts);        

            % [xi_N]:
            J_b_st(:,N_jnts) = xi_R6_body_{N_jnts};
            g_st_SE3_n_i = eye(4); % init
            
            % [xi_N-1' ... xi_1']:
            for i=N_jnts - 1: -1: 1
                % -> compute current jacobian:
                % compute transformation for 1 ... i-1:
                g_st_SE3_n_i = g_st_SE3_n_i * neg_exp_xi_theta_in_SE3_body_{i+1}; 
                % compute adjoint for N ... i:
                Ad_1_i = Lie.Ad_SE3_from_SE3(g_st_SE3_n_i);
                % compute jacobian @ xi_i:
                J_b_st(:,i) =  Ad_1_i * xi_R6_body_{i};
            end
        end
        function J_b_st = compute_Body_Jacobian_from_Spatial(xi_R6_1_, exp_xi_theta_in_SE3_1_, g0_st_SE3_)
            N_jnts = length(xi_R6_1_);
            
            % J = [xi_1, xi_2', ... , xi_N' ]:
            J_b_st = zeros(6,N_jnts);
            g_st_SE3_i_n = eye(4);
            
            % [xi_N' ... xi_1']:
            for i=N_jnts:-1:1
                % -> compute current jacobian:
                g_st_SE3_i_n = exp_xi_theta_in_SE3_1_{i} * g_st_SE3_i_n; 
                Ad_i_n = Lie.inv_Ad_SE3_from_SE3(g_st_SE3_i_n * g0_st_SE3_);
                % compute jacobian @ xi_i:
                J_b_st(:,i) =  Ad_i_n * xi_R6_1_{i};
            end
        end
        %% Shortcuts:
        function g_st_SE3 = compute_Spatial_FK_from_theta(xi_R6_1_, theta_, g0_st_SE3)
            exp_xi_theta_in_SE3_1_ = OpenChain.batch_screw_SE3_from_twist_angle_R6xR(xi_R6_1_, theta_);
            g_st_SE3 = OpenChain.compute_Spatial_FK(exp_xi_theta_in_SE3_1_, g0_st_SE3);
        end
        %% OLD
%         function [G_SE3_b_st_, J_6x6_b_st_] = compute_Body_from_SE3(G_SE3_tool_frame, xi_R6_0_, exp_xi_theta_in_SE3_)
%             % init:
%             N_jnts = length(exp_xi_theta_in_SE3_); 
%             J_6x6_b_st_ = []; % joints starting from N ... 1 
%             G_SE3_b_st_ = {}; % frames starting from N ... 0
%             
%             %% Compute Body (Backward):
%             for i=1:N_jnts
%                 % -> cumulate N's transformation:
%                 G_SE3_b_st_{i+1} = G_SE3_b_st_{i} * exp_xi_theta_in_SE3_{N_jnts-i+1};
%                 % -> compute N's jacobian:
%                 J_6x6_b_st_{i} = Lie.inv_Ad_SE3_from_SE3(G_SE3_b_st_{i}) * xi_R6_0_{N_jnts-i+1};
% 
%                 % -> compute current jacobian:
%                 Ad_0_i = Lie.Ad_SE3_from_SE3(G_SE3_s_st_{i}); % frame i-1:0
%                 J_6x6_s_st_ =  Ad_0_i * xi_R6_0_{i};          % joint   i:1
% 
%                 % -> cumulate current transformation:
%                 G_SE3_s_st_{i+1} = G_SE3_s_st_{i} * exp_xi_theta_in_SE3_{i}; % frame i:1
%             end
%         end
    end
end
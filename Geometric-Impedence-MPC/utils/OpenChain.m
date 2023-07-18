classdef OpenChain
    methods(Static)
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
    end
end
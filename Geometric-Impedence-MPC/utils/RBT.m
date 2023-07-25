classdef RBT
    methods(Static)
    % [SE3]:
        function mat_inv_SE3 = inverse_SE3(mat_SE3)
            mat_inv_SE3 = Lie.inverse_SE3(mat_SE3);
        end
        function mat_SE3 = SE3_from_SO3xR3(R_SO3, p_R3)
            mat_SE3 = Lie.SE3_from_SO3xR3(R_SO3, p_R3);
        end
    % [SO3 x R3]:
        function mat_inv_SE3 = inverse_SO3xR3(R_SO3, p_R3)
            mat_inv_SE3 = Lie.inverse_SO3xR3(R_SO3, p_R3);
        end
        function [R_SO3, p_R3] = SO3xR3_from_SE3(mat_SE3)
            [R_SO3, p_R3] = Lie.SO3xR3_from_SE3(mat_SE3);
        end
    % [Twist Coordinate]:
        function xi = twist_R6_from_axis_point(w_R3, q_R3)
            % input: 
            %   - w \in R^3, axis of rotation, unit vector
            validate.if_unitVector("w_R3", w_R3);
            %   - q \in R^3, point on axis of rotation
            validate.if_dimension("q_R3", q_R3, [3,1]);
            
            % return: xi \in se(3) coordinates vector R^{6}
            %   - rotation about a line
            v_R3 = - Lie.hat_so3_from_R3(w_R3) * q_R3; % (Murray 2.26)
            xi = [ v_R3; w_R3 ]; % twist coordinates
            % hat_xi = Lie.hat_se3_from_R6(xi); % twist
            % Usage in lie algebra:
            % [v; 0] = [dP/dt; 0] = HAT{xi} [P; 1] 
        end
        function xi = twist_R6_from_axis_point_pitch(w_R3, q_R3, h_R)
            % input: 
            %   - w \in R^3, axis of rotation, unit vector
            validate.if_unitVector("w_R3", w_R3);
            %   - q \in R^3, point on axis of rotation
            validate.if_dimension("q_R3", q_R3, [3,1]);
            
            % return: xi \in se(3) coordinates vector R^{6}
            %   - rotation about a line + translation along the line
            v_R3 = - Lie.hat_so3_from_R3(w_R3) * q_R3 + h_R * w_R3; % (Murray 2.26)
            xi = [ v_R3; w_R3 ]; % twist coordinates
        end
    % [Screw Motion]:
        % equivalent:
            % 1. screw_SE3_from_xi_theta
            % 2. screw_SE3_from_axis_point_theta (pure rotation only)
            % 3. screw_SE3_from_axis_point_theta_pitch
        function mat_SE3 = screw_SE3_from_xi_theta(xi_vw, t_R)
            % input: xi_vw \in R6, t_R \in R
            % compute screw motion, pure rotation
            validate.if_dimension("xi_R6", xi_vw, [6,1]);
            mat_SE3 = Lie.exp_map_SE3_from_R6xR(xi_vw, t_R);
            % return: mat_SE3 \in SE(3) \subset R^{4x4} 
        end
        function mat_SE3 = screw_SE3_from_axis_point_theta(w_R3, q_R3, t_R)
            % Pure Rotation, Revolute Joint.
            mat_SE3 = RBT.screw_SE3_from_axis_point_theta_pitch(w_R3, q_R3, t_R, 0);
        end
        function mat_SE3 = screw_SE3_from_axis_point_theta_pitch(w_R3, q_R3, t_R, h_R)
            % input: 
            %   - w \in R^3, axis of rotation, unit vector
            validate.if_unitVector("w_R3", w_R3);
            validate.if_dimension("w_R3", w_R3, [3,1]);
            %   - q \in R^3, point on axis of rotation
            validate.if_dimension("q_R3", q_R3, [3,1]);
            % (Murray 2.40) <--- reduced from exp_map (2.36)
            R = Lie.rodrigues_SO3_from_unit_R3xR(w_R3,t_R);
            if h_R == 0
                t = (eye(3)-R)*q_R3;
            else
                t = (eye(3)-R)*q_R3 + h_R*w_R3*t_R; % translation vector
            end
            mat_SE3 = [R, t; zeros(1,3), 1];
            % return: mat_SE3 \in SE(3) \subset R^{4x4} 
        end
        function mat_SE3 = screw_SE3_pure_translation(t_R, v_R3)
            validate.if_unitVector("v_R3", v_R3);
            validate.if_dimension("v_R3", v_R3, [3,1]);
            mat_SE3 = [
                eye(3), t_R*v_R3;
                zeros(1,3), 1
            ];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Modern Robotics:
        function R_SO3 = project_to_SO3_from_nearSO3(mat_near_SO3)
            [U, S, V] = svd(mat_near_SO3);
            R = U * V';
            if det(R) < 0
                % In this case the result may be far from mat.
                R_SO3 = [R(:, 1: 2), -R(:, 3)];
            end
        end
        function mat_SE3 = project_to_SE3_from_nearSE3(mat_near_SE3)
            mat_SE3 = [ RBT.project_to_SO3_from_nearSO3(mat_near_SE3(1: 3, 1: 3)), mat_near_SE3(1: 3, 4); zeros(1, 3), 1];
        end
    end
end
classdef Lie
    properties(Constant)
        INFINITY = 1e+9;
    end
    methods(Static)
    % VEE and HAT operators for Lie algebra:
        function vec_vw = vee_R6_from_se3(mat_R4x4)
            % input: mat \in R^{4x4} 
            validate.if_dimension("mat_R4x4", mat_R4x4, [4,4]);
            % v:
            vec_vw(1:3) = mat_R4x4(1:3,4);
            % w:
            vec_vw(4) = mat_R4x4(3,2);
            vec_vw(5) = mat_R4x4(1,3);
            vec_vw(6) = mat_R4x4(2,1);
            % return: vec \in R^{6}
        end
        function vec = vee_R3_from_so3(mat_R3x3)
            % input: mat \in R^{3x3} 
            validate.if_dimension("mat_R3x3", mat_R3x3, [3,3]);
            vec(1) = mat_R3x3(3,2);
            vec(2) = mat_R3x3(1,3);
            vec(3) = mat_R3x3(2,1);
            % return: vec \in R^{3}
        end
        function mat_R4x4 = hat_se3_from_R6(vec_vw)
            % input: vec \in se(3) coordinates vector R^{6}
            validate.if_dimension("(v,w) \in se(3) coordinates", vec_vw, [6,1]);
            mat_R4x4 = [
                Lie.hat_so3_from_R3(vec_vw(4:6)), vec_vw(1:3); ...
                zeros(1,4)
            ];
            % return: mat \in se(3) \subset R^{4x4} 
        end
        function mat_R3x3 = hat_so3_from_R3(vec)
            % input: vec \in R^{3}
            validate.if_dimension("vec", vec, [3,1]);
            w = vec;
            mat_R3x3 = [
                        0, -w(3),  w(2); ...
                     w(3),     0, -w(1); ...
                    -w(2),  w(1),    0
            ];
            % return: mat \in R^{3x3} 
        end
    % [SE3]:
        function mat_inv_SE3 = inverse_SE3(mat_SE3)
            [R_SO3, p_R3] = Lie.SO3xR3_from_SE3(mat_SE3);
            mat_inv_SE3 = Lie.inverse_SO3xR3(R_SO3, p_R3);
        end
        function mat_SE3 = SE3_from_SO3xR3(R_SO3, p_R3)
            mat_SE3 = [
                R_SO3, p_R3;
                zeros(1, 3), 1
            ];
        end
    % [SO3 x R3]:
        function mat_inv_SE3 = inverse_SO3xR3(R_SO3, p_R3)
            mat_inv_SE3 = [
                R_SO3.', -R_SO3.' * p_R3;
                zeros(1, 3), 1
            ];
        end
        function [R_SO3, p_R3] = SO3xR3_from_SE3(mat_SE3)
            R_SO3 = mat_SE3(1:3, 1:3);
            p_R3 = mat_SE3(1:3, 4);
        end
    % [Rodrigue]
        function exp_SO3 = rodrigues_SO3_from_R3(w_R3)
            % > input: w_R3 \in so(3) coordinate \subset R^{3} vector
            % - [rodrigues' formula]
            % - so3 ---> SO3
            % < output: mat \in SO(3) rotation matrix R^{3x3}
            validate.if_dimension("w_R3", w_R3, [3,1]);
            % pre-compute const, (minimize division operations)
            w_abs = norm(w_R3,2);
            
            exp_SO3 = eye(3);
            if w_abs ~= 0
                inv_w_abs = 1 / w_abs;
                inv_w_abs_sqr = inv_w_abs * inv_w_abs;
                hat_w_so3 = Lie.hat_so3_from_R3(w_R3);
                
                % [ rodrigues' formula std ]:
                exp_SO3 = exp_SO3 + sin(w_abs) * hat_w_so3 * inv_w_abs + (1 - cos(w_abs)) * hat_w_so3^2 * inv_w_abs_sqr;
            end
            
            % [ validation ]:
            validate.if_SO("rodrigues_SO3_from_R3", exp_SO3);
        end
        function exp_SO3 = rodrigues_SO3_from_so3(hat_w_so3)
            % > input: hat_w_so3 \in so(3) \subset R^{3x3} skew-symmetric matrix
            % - [rodrigues' formula]
            % - so3 ---> SO3
            % < output: mat \in SO(3) rotation matrix R^{3x3}
            validate.if_dimension("hat_w_so3", hat_w_so3, [3,3]);
            validate.if_so_("hat_w_so3", hat_w_so3);
            
            % pre-compute const, (minimize division operations)
            w_R3 = Lie.vee_R3_from_so3(hat_w_so3);
            w_abs = norm(w_R3,2);
            inv_w_abs = 1 / w_abs;
            inv_w_abs_sqr = inv_w_abs * inv_w_abs;
            
            % [ rodrigues' formula std ]:
            exp_SO3 = eye(3) + sin(w_abs) * hat_w_so3 * inv_w_abs + (1 - cos(w_abs)) * hat_w_so3^2 * inv_w_abs_sqr;
            
            % [ validation ]:
            validate.if_SO("rodrigues_SO3_from_so3", exp_SO3);
        end
        function exp_SO3 = rodrigues_SO3_from_unit_R3xR(w_R3, t_R)
            % > input: w_R3 \in so(3) coordinate \subset R^{3} unit vector, t_R \in R scalar 
            % - [rodrigues' formula]: axis-angle representation of rotation matrix
            % - so3 ---> SO3
            % < output: mat \in SO(3) rotation matrix R^{3x3}
            validate.if_dimension("w_R3", w_R3, [3,1]);
            validate.if_unitVector("w_R3", w_R3);

            % pre-compute const, (minimize division operations)
            exp_SO3 = eye(3);
            hat_w_so3 = Lie.hat_so3_from_R3(w_R3);
            
            % [ rodrigues' formula unit vector x scalar ]:
            exp_SO3 = exp_SO3 + sin(t_R) * hat_w_so3 + (1 - cos(t_R)) * hat_w_so3^2;
            
            % [ validation ]:
            validate.if_SO("rodrigues_SO3_from_unit_R3xR", exp_SO3);
        end 
        function mat_SE3 = exp_map_SE3_from_R6(xi_vw)
            % NOTE: may result different from exp_map_SE3_from_R6xR, but this is general form:
            % twist coordinates vector \in R^{6} (v,w)
            validate.if_dimension("xi_vw", xi_vw, [6,1]);
            v_R3 = xi_vw(1:3);
            w_R3 = xi_vw(4:6);
            
            % pre-compute const, (minimize division operations)
            w_abs = norm(w_R3);
            
            if w_abs ~= 0
                inv_w_abs = 1 / w_abs;
                inv_w_abs_sqr = inv_w_abs * inv_w_abs;
                hat_w_so3 = Lie.hat_so3_from_R3(w_R3);
                
                % [ rodrigues' formula std ]:
                R = eye(3) + sin(w_abs) * hat_w_so3 * inv_w_abs + (1 - cos(w_abs)) * hat_w_so3^2 * inv_w_abs_sqr;
        
                % [ validation ]:
                % validate.if_SO("{R | mat_SE3(R,p)}", R);
            
                % [ transplation ]:
                p = inv_w_abs_sqr * ((eye(3) - R) * hat_w_so3 + w_R3 * w_R3') * v_R3;
                
                % [ OUTPUT ]:
                mat_SE3 = [R, p; zeros(1,3), 1];
            else
                mat_SE3 = [eye(3), v_R3; zeros(1,3), 1];
            end
        end
        function exp_SE3 = exp_map_SE3_from_R6xR(xi_vw, t_R)
            % NOTE: may result different from exp_map_SE3_from_R6,
            %       but this is a more efficient way with less numerical error from small angles.
            % This is a special case of exp_map_SE3_from_R6, where w_R3 is unit vector.
            %
            % input: xi_vw: twist coordinates vector \in R^{6} (v,w)
            validate.if_dimension("xi_vw", xi_vw, [6,1]);
            v_R3 = xi_vw(1:3);
            w_R3 = xi_vw(4:6);
            
            validate.if_unitVector("w_R3", w_R3);
            
            % pre-compute const, (minimize division operations):
            if norm(w_R3) == 0
                exp_SE3 = [eye(3), v_R3 * t_R; zeros(1,3), 1]; % (Murray 2.32)
            else
                % [ rodrigues' formula std ]:
                R = Lie.rodrigues_SO3_from_unit_R3xR(w_R3, t_R);
                
                % [ transplation ]:
                hat_w_so3 = Lie.hat_so3_from_R3(w_R3);
                % validate.if_inRangeStrict("t_R", t_R, (0, 2*pi)); 
                % --> A_mat(.) non-singular for all (0,2*pi)
                A_mat = ((eye(3) - R) * hat_w_so3 + w_R3 * w_R3.' * t_R); % [Murray 2.38]
                p = A_mat * v_R3; % (Murray 2.36)
                
                % [ OUTPUT ]:
                exp_SE3 = [R, p; zeros(1,3), 1];
            end
        end
    % [Log]
        function mat_so3 = log_so3_from_SO3(R_SO3)
            validate.if_SO("R_SO3", R_SO3);
            acosinput = (trace(R_SO3) - 1) / 2;
            if acosinput >= 1
                mat_so3 = zeros(3);
            elseif acosinput <= -1
                if ~ (norm(1 + R(3, 3)) < 1e-6)
                    omg = (1 / sqrt(2 * (1 + R(3, 3)))) ...
                        * [R(1, 3); R(2, 3); 1 + R(3, 3)];
                elseif ~ (norm(1 + R(2, 2)) < 1e-6)
                    omg = (1 / sqrt(2 * (1 + R(2, 2)))) ...
                        * [R(1, 2); 1 + R(2, 2); R(3, 2)];
                else
                    omg = (1 / sqrt(2 * (1 + R(1, 1)))) ...
                        * [1 + R(1, 1); R(2, 1); R(3, 1)];
                end
                mat_so3 = Lie.hat_so3_from_R3(pi * omg);
            else
                theta = acos(acosinput);
                mat_so3 = theta * (1 / (2 * sin(theta))) * (R_SO3 - R_SO3');
            end
        end
        function mat_se3 = log_se3_from_SE3(mat_SE3)
            % MatrixLog6
            [R_SO3, p_R3] = Lie.SO3xR3_from_SE3(mat_SE3);
            hat_omega = Lie.log_so3_from_SO3(R_SO3);
            if isequal(hat_omega, zeros(3))
                mat_se3 = [zeros(3), T(1: 3, 4); 0, 0, 0, 0];
            else
                theta = acos((trace(R_SO3) - 1) / 2);
                mat_se3 = [hat_omega, (eye(3) - hat_omega / 2 ...
                                + (1 / theta - cot(theta / 2) / 2) ...
                                    * hat_omega * hat_omega / theta) * p_R3;
                        0, 0, 0, 0];    
            end
        end
    % [Adjoint]
        function mat = Ad_SE3_from_SE3(mat_SE3)
            % input: mat_SE3 \in SE(3) \subset R^{4x4} 
            R_SO3 = mat_SE3(1:3,1:3);
            p_R3 = mat_SE3(1:3,4);
            validate.if_SO("{R | mat_SE3(R,p)}", R_SO3);
            validate.if_dimension("mat_SE3", mat_SE3, [4,4]);
            hat_p = Lie.hat_so3_from_R3(p_R3);
            mat = [
                R_SO3(:,:), hat_p * R_SO3(:,:); ...
                zeros(3), R_SO3(:,:)
            ];
            % return: mat \in R^{6x6}
        end
        function mat = inv_Ad_SE3_from_SE3(mat_SE3)
            % input: mat_SE3 \in SE(3) \subset R^{4x4} 
            R_SO3 = mat_SE3(1:3,1:3);
            p_R3 = mat_SE3(1:3,4);
            validate.if_SO("{R | mat_SE3(R,p)}", R_SO3);
            validate.if_dimension("mat_SE3", mat_SE3, [4,4]);
            hat_p = Lie.hat_so3_from_R3(p_R3);
            mat = [
                R_SO3(:,:)', -R_SO3(:,:)' * hat_p; ...
                zeros(3), R_SO3(:,:)'
            ]; % [Murray 2.58]
            % return: mat \in R^{6x6}
        end
        function mat = Ad_SE3_from_SO3xR3(R_SO3, p_R3)
            % input: mat_SE3 \in SE(3) \subset R^{4x4} 
            validate.if_SO("{R | mat_SE3(R,p)}", R_SO3);
            validate.if_dimension("R_SO3", R_SO3, [3,3]);
            validate.if_dimension("p_R3", p_R3, [3,1]);
            hat_p = Lie.hat_so3_from_R3(p_R3);
            mat = [
                R_SO3(:,:), hat_p * R_SO3(:,:); ...
                zeros(3), R_SO3(:,:)
            ]; % [Murray 2.58]
            % return: mat \in R^{6x6}
        end
        function mat = inv_Ad_SE3_from_SO3xR3(R_SO3, p_R3)
            % input: mat_SE3 \in SE(3) \subset R^{4x4} 
            validate.if_SO("{R | mat_SE3(R,p)}", R_SO3);
            validate.if_dimension("R_SO3", R_SO3, [3,3]);
            validate.if_dimension("p_R3", p_R3, [3,1]);
            hat_p = Lie.hat_so3_from_R3(p_R3);
            mat = [
                R_SO3(:,:)', -R_SO3(:,:)' * hat_p; ...
                zeros(3), R_SO3(:,:)'
            ];
            % return: mat \in R^{6x6}
        end
        function mat = ad_se3_from_R6(vec_vw)
            % input: vec_vw \in se(3) coordinates vector R^{6}
            validate.if_dimension("(v,w)", vec_vw, [6,1]);
            v = vec_vw(1:3,:);
            w = vec_vw(4:6,:);
            hat_w = Lie.hat_so3_from_R3(w);
            hat_v = Lie.hat_so3_from_R3(v);
            mat = [
                hat_w(:,:), hat_v(:,:); ...
                zeros(3), hat_w(:,:)
            ];
        end
    % [SE3 products]
        function mat = prodLeft_SE3_from_SE3xk(mats)
            mat = mats{1};
            % N x N-1 x ... x 1
            for i = 2:length(mats)
                mat = mat * mats{i};
            end
        end
        function mat = prodRight_SE3_from_SE3xk(mats)
            mat = mats{1};
            % 1 x ... x N
            for i = 2:length(mats)
                mat = mats{i} * mat;
            end
        end
    % [SE3 ops]
        function d = distance_to_SO3(mat_SO3)
            % *** CHAPTER 3: RIGID-BODY MOTIONS ***        
            if det(mat_SO3) > 0
                d = norm(mat_SO3' * mat_SO3 - eye(3), 'fro');
            else
                d = Lie.INFINITY; % infty
            end
        end
        function d = distance_to_SE3(mat_SE3)
            % *** CHAPTER 3: RIGID-BODY MOTIONS ***
            R_SO3 = mat_SE3(1:3,1:3);
            if det(R_SO3) > 0
                d = norm([R' * R, [0; 0; 0]; mat_SE3(4, :)] - eye(4), 'fro');
            else
                d = Lie.INFINITY; % infty
            end
        end
    % [Lie Bracket]
        function vec = lie_bracket_from_R6xR6(xi_1, xi_2) 
            % strictly [xi_1, xi_2] = -[xi_2, xi_1]
            xi_1_hat = Lie.hat_se3_from_R6(xi_1);
            xi_2_hat = Lie.hat_se3_from_R6(xi_2);
            vec = Lie.lie_bracket_from_se3xse3(xi_1_hat, xi_2_hat);
        end
        function vec = lie_bracket_from_se3xse3(xi_1_hat, xi_2_hat) 
            vec = Lie.vee_R6_from_se3(xi_1_hat * xi_2_hat - xi_2_hat * xi_1_hat);
        end
    end
end
function mat_SE3 = screw_SE3_from_xi_theta(xi_R6, t_R)
    % input: xi_R6 \in R6, t_R \in R
    % compute screw motion, pure rotation
    validate.if_dimension("xi_R6", xi_R6, [6,1]);
    
    q = xi_R6(1:3);
    R = rodrigues_SO3_from_R3xR(xi_R6(4:6),t_R);
    % reshape to [:,:,i]
    mat_SE3 = [R, (eye(3)-R)*q; zeros(1,3), 1]; % (Murray 2.40) <--- reduced from exp_map (2.36)
    % return: mat_SE3 \in SE(3) \subset R^{4x4} 
end
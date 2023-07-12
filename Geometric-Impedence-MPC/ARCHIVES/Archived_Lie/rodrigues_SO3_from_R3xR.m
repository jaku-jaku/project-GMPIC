function exp_mat = rodrigues_SO3_from_R3xR(w_R3, t_R)
    % > input: w_R3 \in so(3) coordinate \subset R^{3} unit vector, t_R \in R scalar 
    % - [rodrigues' formula]: axis-angle representation of rotation matrix
    % - so3 ---> SO3
    % < output: mat \in SO(3) rotation matrix R^{3x3}
    validate.if_dimension("w_R3", w_R3, [3,1]);
    validate.if_unitVector("w_R3", w_R3);
    % validate.if_inRangeStrict("t_R", t_R, [0, 2*pi]); % --> A(.) non-singular

    % pre-compute const, (minimize division operations)
    exp_mat = eye(3);
    hat_w_so3 = hat_so3(w_R3);
    
    % [ rodrigues' formula unit vector x scalar ]:
    exp_mat = exp_mat + sin(t_R) * hat_w_so3 + (1 - cos(t_R)) * hat_w_so3^2;
    
    % [ validation ]:
    validate.if_SO3("R_rodrigue_SO3", exp_mat);
end
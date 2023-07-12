function mat = hat_se3(vec)
    % input: vec \in se(3) coordinates vector R^{6xk}
    validate.if_dimension("vec \in se(3) coordinates", vec, [6,1]);
    v = vec(1:3);
    w = vec(4:6);
    mat = [
            0, -w(3),  w(2), v(1); ...
         w(3),     0, -w(1), v(2); ...
        -w(2),  w(1),     0, v(3); ...
        0,0,0,0
    ];
    % return: mat \in se(3) \subset R^{4x4xk} 
end
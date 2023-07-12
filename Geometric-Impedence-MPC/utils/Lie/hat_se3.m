function mat = hat_se3(vec)
    % input: vec \in se(3) coordinates vector R^{6xk}
    [ n, k ] = size(vec);
    assert(n == 6, "[ERR] hat_se3: vec \in se(3) coordinates vector R^{6} (v,w)!");
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
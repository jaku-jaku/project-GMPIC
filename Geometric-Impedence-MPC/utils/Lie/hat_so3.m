function mat = hat_so3(vec)
    % input: vec \in so(3) coordinates vector R^{3}
    [ n, k ] = size(vec);
    assert(n == 3, "[ERR] hat_so3: vec must be  R^{3} vectors!");
    w = vec;
    mat = [
                0, -w(3),  w(2); ...
             w(3),     0, -w(1); ...
            -w(2),  w(1),    0
    ];
    % return: mat \in so(3) \subset R^{3x3} 
end
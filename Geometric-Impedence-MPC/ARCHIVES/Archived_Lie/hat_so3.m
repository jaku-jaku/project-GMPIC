function mat = hat_so3(vec)
    % input: vec \in so(3) coordinates vector R^{3}
    validate.if_dimension("vec", vec, [3,1]);
    w = vec;
    mat = [
                0, -w(3),  w(2); ...
             w(3),     0, -w(1); ...
            -w(2),  w(1),    0
    ];
    % return: mat \in so(3) \subset R^{3x3} 
end
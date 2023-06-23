function mat = hat_so3(vec)
    % input: vec \in so(3) coordinates vector R^{3xk}
    [n, m] = size(vec);
    assert(n == 3, "[ERR] hat: vec must be  R^{3xk} vectors!");
    mat=[];
    w = vec;
    for i = 1:m
        mat = [mat,  [
            0, -w(3,i),  w(2,i); ...
            w(3,i), 0, -w(1,i); ...
            -w(2,i), w(1,i), 0
        ]];
    end
    % reshape to [:,:,i]
    mat = reshape(mat, 3, 3, m);
    % return: mat \in so(3) \subset R^{3x3xk} 
end
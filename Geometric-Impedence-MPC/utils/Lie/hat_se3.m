function mat = hat_se3(vec)
    % input: vec \in se(3) coordinates vector R^{6xk}
    [n, m] = size(vec);
    assert(n == 6, "[ERR] ad_se3_: vec \in se(3) coordinates vector R^{6xk} (v,w)!");
    mat=[];
    v = vec(1:3,:);
    w = vec(4:6,:);
    for i = 1:m
        mat = [mat,  [
            0, -w(3,i),  w(2,i), v(1,i); ...
            w(3,i), 0, -w(1,i), v(2,i); ...
            -w(2,i), w(1,i), 0, v(3,i); ...
            0,0,0,0
        ]];
    end
    % reshape to [:,:,i]
    mat = reshape(mat, 4, 4, m);
    % return: mat \in se(3) \subset R^{4x4xk} 
end
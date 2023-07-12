function mat = batch_ad_se3_from_R6(vec)
    % input: vec \in se(3) coordinates vector R^{6xk}
    [n, m] = size(vec);
    assert(n == 6, "[ERR] ad_se3_: vec \in se(3) coordinates vector R^{6xk} (v,w)!");
    v = vec(1:3,:);
    w = vec(4:6,:);
    hw = hat_so3(w);
    hv = hat_so3(v);
    mat=[];
    for i = 1:m
        mat = [ mat, [
            hw(:,:,i), zeros(3); ...
            hv(:,:,i), hw(:,:,i)
        ]];
    end
    mat = reshape(mat, 6, 6, m);
    % return: mat \in R^{6x6xk}
end
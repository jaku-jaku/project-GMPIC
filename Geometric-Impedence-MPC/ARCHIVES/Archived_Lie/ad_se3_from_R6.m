function mat = ad_se3_from_R6(vec)
    % input: vec \in se(3) coordinates vector R^{6}
    n = size(vec);
    assert(n == 6, "[ERR] ad_se3_: vec \in se(3) coordinates vector R^{6} (v,w)!");
    v = vec(1:3,:);
    w = vec(4:6,:);
    hw = hat_so3(w);
    hv = hat_so3(v);
    mat = [
        hw(:,:), zeros(3); ...
        hv(:,:), hw(:,:)
    ];
end
function mat = Ad_SE3_from_SE3_or_SO3xR3(arg1,arg2)
    if nargin == 1
        % input: mat_SE3 \in SE(3) \subset R^{4x4} 
        mat_SE3 = arg1;
        R = mat_SE3(1:3,1:3);
        p = mat_SE3(1:3,4);
        [n, m] = size(mat_SE3);
        assert(prod([n,m] == [4,4]), "[ERR] Ad_SE3: mat_SE3 \in SE(3) \subset R^{4x4}!");
    else
        % input: R \in SO(3), p \in \R^{3} 
        R=arg1;
        p=arg2;
        [n, m] = size(R);
        [n2] = size(p);
        assert(prod([n,m,n2] == [3,3,3]), "[ERR] R \in SO(3), p \in \R^{3} !");
    end
    
    mat = [
        R(:,:), zeros(3); ...
        hat_so3(p) * R(:,:), R(:,:)
    ];
    % return: mat \in R^{6x6xk}
end
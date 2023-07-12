function mat = Ad_SE3_from_SE3_or_SO3xR3(arg1,arg2)
    if nargin == 1
        % input: mat_SE3 \in SE(3) \subset R^{4x4xk} 
        mat_SE3 = arg1;
        R = mat_SE3(1:3,1:3,:);
        p = mat_SE3(1:3,4,:);
        [n, m, k] = size(mat_SE3);
        assert(prod([n,m] == [4,4]), "[ERR] Ad_SE3: mat_SE3 \in SE(3) \subset R^{4x4xk}!");
    else
        % input: R \in SO(3)xk, p \in \R^{3xk}  
        R=arg1;
        p=arg2;
        [n, m, k] = size(R);
        [n, k2] = size(p);
        assert(prod([n,m] == [3,3]), "[ERR] R \in SO(3)xk, p \in \R^{3xk} !");
        assert(k==k2, "R, p should have same length!");
    end
    
    hp = hat_so3(p);
    mat=[];
    
    for i = 1:k
        mat = [ mat, [
            R(:,:,i), zeros(3); ...
            hp(:,:,i) * R(:,:,i), R(:,:,i)
        ]];
    end
    mat = reshape(mat, 6, 6, k);
    % return: mat \in R^{6x6xk}
end
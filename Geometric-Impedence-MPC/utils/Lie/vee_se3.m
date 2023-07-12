function vec = vee_se3(mat)
    % input: mat \in se(3) \subset R^{4x4} 
    [n, m] = size(mat);
    assert(prod([n,m] == [4,4]), "[ERR] vee: mat must be R^{4x4} matrices!");
    % v:
    vec(1:3) = mat(1:3,4);
    % w:
    vec(4) = -mat(2,3);
    vec(5) =  mat(1,3);
    vec(6) = -mat(1,2);
    % return: vec \in se(3) coordinates vector R^{6}
end

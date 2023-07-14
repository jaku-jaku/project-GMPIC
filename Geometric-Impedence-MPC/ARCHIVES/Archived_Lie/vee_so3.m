function vec = vee_so3(mat)
    % input: mat \in so(3) \subset R^{3x3} 
    [n, m] = size(mat);
    assert(prod([n,m] == [3,3]), "[ERR] vee: mat must be R^{3x3} matrices!");
    vec(1) = -mat(2,3);
    vec(2) =  mat(1,3);
    vec(3) = -mat(1,2);
    % return: vec \in so(3) coordinates vector R^{3}
end

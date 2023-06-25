function mat_SE3 = exp_map_SE3(xi)
    [n, m] = size(xi);
    assert(n == 6, "[ERR] ad_se3_: xi \in se(3) coordinates vector R^{6xk} (v,w)!");
    v = xi(1:3,:);
    w = xi(4:6,:);
    
    % Optimization for numerical only:
        % ZERO_TOLERANCE_THRESHOLD = 1e-6;
        % For w = 0
        % if norm(w - 0, 1) < ZERO_TOLERANCE_THRESHOLD
        %     % [ efficient closed form expression ]: \cite{parkGeometricIntegrationEuclidean2005}
        %     mat_SE3 = [eye(3), v; 0, 0, 0, 1];
        %     % else
        % else
    
    % General:
    % [ efficient closed form expression ]: \cite{parkGeometricIntegrationEuclidean2005}
    Omega = rodrigues_SO3_from_R3(w);
    mat_12 = 1 / norm(w, 2) * ((eye(3) - Omega) * cross(w, v) + w * w' * v);
    mat_SE3 = [Omega, mat_12; zeros(1,3), 1];
end
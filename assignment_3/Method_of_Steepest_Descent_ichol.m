function [x,niters] = Method_of_Steepest_Descent_ichol(A,b)
    L = sparse( ichol(sparse(A), struct('type','ict','droptol',1e-2,'michol','off')));
    x = ones(size(b, 1), 1);
    r = b - A * x;
    niters = 0;
    stop_cond = sqrt(eps) * norm(b);
    while norm(r) >= stop_cond
        z = L \ r;
        p = L' \ z;
        q = A * p;
        alpha = dot(p, r) / dot(p, q);
        x = x + alpha * p;
        r = r - alpha * q;
        niters = niters + 1;
    end
end
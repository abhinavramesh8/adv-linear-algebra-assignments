function [x,niters] = PCG(A,b)
    L = sparse( ichol(sparse(A), struct('type','ict','droptol',1e-2,'michol','off')));
    x = ones(size(b, 1), 1);
    r = b - A*x;
    niters = 0;
    stop_cond = sqrt(eps) * norm(b);
    while norm(r) >= stop_cond
        t = L \ r;
        z = L' \ t;
        rz = dot(r, z);
        if niters == 0
            p = z;
        else
            gamma = rz / oldrz;
            p = z + gamma * p;
        end
        q = A * p;
        alpha = rz / dot(p, q);
        x = x + alpha * p;
        r = r - alpha * q;
        oldrz = rz;
        niters = niters + 1;
    end
end
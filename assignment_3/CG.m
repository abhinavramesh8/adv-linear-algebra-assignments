function [x,niters] = CG(A,b)
    x = ones(size(b, 1),1);
    r = b - A*x;
    niters = 0;
    stop_cond = sqrt(eps) * norm(b);
    while norm(r) >= stop_cond
        rsq = dot(r,r);
        if niters == 0
            p = r;
        else
            gamma = rsq / oldrsq;
            p = r + gamma*p;
        end
        ap = A*p;
        ptap = dot(p, ap);
        alpha = rsq / ptap;
        x = x + alpha*p;
        r = r - alpha*ap;
        oldrsq = rsq;
        niters = niters + 1;
    end
end
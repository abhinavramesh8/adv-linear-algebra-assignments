function [x,niters] = Method_of_Steepest_Descent(A,b)
    x = ones(size(b, 1), 1);
    niters = 0;
    r = b - A*x;
    stop_cond = sqrt(eps) * norm(b);
    while norm(r) >= stop_cond
        p = r;
        q = A*p;
        alpha = dot(p,r) / dot(p,q);
        x = x + alpha*p;
        r = r - alpha*q;
        niters = niters + 1;
    end
end
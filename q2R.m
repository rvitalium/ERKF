function R = q2R(q)
    q = reshape(q,[4,1]);
    w = q(1);
    v = q(2:4);
    
    R = (w^2 - v'*v)*eye(3) + 2*v*v' + 2*w*skew(v);
end
function q_conj = quat_conj(q)
    s = size(q);
    q = reshape(q,1,[]);
    q_conj = reshape([q(1),-q(2:4)],s);
end
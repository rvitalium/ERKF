% error state kinematics

function P_new = errStateKin(x,u,dt,P,Qi,xind,uind)
    
    q = vec2quat(u(uind.w)-x(xind.wb),dt);
    R = q2R(q');
    
    %         dth        dwb
    Fx = [      R', -eye(3)*dt;  % dth
           zeros(3),    eye(3)]; % dwb
    
    Fi = blkdiag(eye(3),eye(3));
    
    P_new = Fx*P*Fx' + Fi*Qi*Fi';
    
end
% nominal state kinematics

function x_new = nomStateKin(x,u,dt,xind,uind)
    
    % orientation
    phi = norm(u(uind.w)-x(xind.wb))*dt;
    U = (u(uind.w)-x(xind.wb))./norm(u(uind.w)-x(xind.wb));
    q = quatMult(x(xind.q),[cos(phi/2);sin(phi/2).*U]);
    q = q./norm(q);
    
    % angular rate gyroscope bias
    wb = x(xind.wb);
    
    % collect new states
    x_new = [q; wb];
    
end
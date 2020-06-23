% update nominal state kinematics

function x_new = nomStateUpd(x,dx,xind,dxind)

    quest = [1;0.5.*dx(dxind.dth)];
    quest = quest./norm(quest);

    x_new = [quatMult(x(xind.q),quest);
             x(xind.wb) + dx(dxind.dwb)];
         
    x_new(xind.q) = x_new(xind.q)./norm(x_new(xind.q));
         
end
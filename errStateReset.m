% update nominal state kinematics

function [dx_new, P_new] = errStateReset(dx,P)

    dx_new = [dx(1:3) - dx(1:3);
              dx(4:6) - dx(4:6);
              dx(7:9) - dx(7:9);
              dx(10:12) - dx(10:12);
              dx(13:15) - dx(13:15)];
    
    G =  blkdiag(eye(6),eye(3)-skew(0.5.*dx(7:9)),eye(6));
    P_new = G*P*G';
    
end
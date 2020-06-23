% update error state using measurements

function [dx_new,P_new] = errStateUpd(x,u,R_,P,g,xind,uind,K)

    q = x(xind.q,end);
    w = q(1);
    v = q(2:4);
    R = q2R(q');
    G = -g./norm(g);
    s = K.static;
    m = K.magnet;
    
    if size(u,2) ~= 1
        if s == 0 && m ~= 0
            mm(:,1) = q2R(x(xind.q,1))*(u(uind.m,1)./norm(u(uind.m,1)));
            mm(:,2) = R*(u(uind.m,end)./norm(u(uind.m,end)));
            psi1 = atan2(mm(2,1),mm(1,1));
            psi2 = atan2(mm(2,2),mm(1,2));
            psi_d = (psi2-psi1)*128;
            if abs(psi_d) < 25*(pi/180)
                m = 0;
                M = u(uind.m,end-1)./norm(u(uind.m,1));
                RM = q2R(x(xind.q,1)');
                M = RM*M;
                zm = u(uind.m,end)./norm(u(uind.m,end));
            end
        end
    end
    
    if ~exist('M','var')
        M = [1;0;0];
        if size(u,2) == 2
            zm = u(uind.m,2)./norm(u(uind.m,2));
        else
            zm = u(uind.m,1)./norm(u(uind.m,1));
        end
        Mk = R*zm;
        Mk(3) = 0;
        Mk = Mk./norm(Mk);
        zm = R'*Mk;
    end
    
    u = u(:,end);
    
    A = char(65+double([logical(s),m]));
    
    switch A
        case {'AA'} % accelerometer and magnetometer data update
            % accelerometer measurement
            za = u(uind.a)./norm(u(uind.a));
            h_xa = R'*G;
            h1 = 2*[w.*G + cross(G,v), v'*G*eye(3) + v*G' - G*v' + w.*skew(G)];

            % magnetometer measurement
            h_xm = R'*M;
            h2 = 2*[w.*M + cross(M,v), v'*M*eye(3) + v*M' - M*v' + w.*skew(M)];

            % combine measurements
            z = [za;zm];
            h_x = [h_xa;h_xm];
            Hx = [h1; h2];
            Hx = [Hx, zeros(6,3)];
            Q_dth = 0.5*[-q(2), -q(3), -q(4);
                          q(1), -q(4),  q(3);
                          q(4),  q(1), -q(2);
                         -q(3),  q(2),  q(1)];
            X_dth = blkdiag(Q_dth,eye(3));
            H = Hx*X_dth;
            
        case {'AB'} % accelerometer data update
            
            z = u(uind.a)./norm(u(uind.a));
            h_x = R'*G;
            h = 2*[w.*G + cross(G,v), v'*G*eye(3) + v*G' - G*v' + w.*skew(G)];
            Hx = [h, zeros(3)];
            Q_dth = 0.5*[-q(2), -q(3), -q(4);
                          q(1), -q(4),  q(3);
                          q(4),  q(1), -q(2);
                         -q(3),  q(2),  q(1)];
            X_dth = blkdiag(Q_dth,eye(3));
            H = Hx*X_dth;
            R_ = R_(1:3,1:3);
            
        case {'BA'} % magnetometer data update
            
            z = zm;
            h_x = R'*M;
            h = 2*[w.*M + cross(M,v), v'*M*eye(3) + v*M' - M*v' + w.*skew(M)];
            Hx = [h, zeros(3)];
            Q_dth = 0.5*[-q(2), -q(3), -q(4);
                          q(1), -q(4),  q(3);
                          q(4),  q(1), -q(2);
                         -q(3),  q(2),  q(1)];
            X_dth = blkdiag(Q_dth,eye(3));
            H = Hx*X_dth;
            R_ = R_(4:6,4:6);
                     
        otherwise % none
            dx_new = zeros(length(x)-1,1);
            P_new = P;
            return;
    end

    % innovation
    I = z - h_x;
    
    % innovation covariance
    S = R_ + H*P*H';
    
    % Kalman gain
    K = P*H'/S;
    
    % update
    dx_new = K*I;
    P_new = (eye(6) - K*H)*P;

end
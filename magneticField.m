function res = magneticField(val,mu,sigma,mu_mag,sigma_mag,x0,x,u,xind,uind,dt)

%{

possible values for res

0 -> magnetically clean environment
1 -> magnetic interference present

%}

alpha_ = 0.05;
m_mag = sqrt(sum(u(uind.m,:).^2,1))';
maDist2 = (mu_mag-m_mag)'/(sigma_mag^2.*eye(length(m_mag)))*(mu_mag-m_mag);
r = maDist2 < chi2inv(1-alpha_,size(maDist2,2));

if r
    if val == 0 % if the sensor is static -> Sabatini's method
        a = u(uind.a,end);
        m = u(uind.m,end);

        for i = 1:size(a,2)
            theta(i,1) = acos((a(:,i)'*m(:,i))/(norm(a(:,i))*norm(m(:,i))));
        end

        maDist2 = (mu-theta)'/(sigma^2.*eye(length(theta)))*(mu-theta);
        r = maDist2 < chi2inv(1-alpha_,size(maDist2,2));

        if r == 1
            res = 0;
        else
            res = 1;
        end

    else % if the sensor is not static -> Kwanmuang's method
        w = u(uind.w,2:end);
        m = u(uind.m,:);

        for i = 1:size(w,2)
            q0 = x0(xind.q,i);
            R0 = q2R(q0');
            q1 = x(xind.q,i);
            R1 = q2R(q1');

            e = quat2eul(q1','ZYZ');
            theta = e(2);
            psi = e(3);

            if theta == 0 || theta == pi
                res = 1;
                return;
            else
                psi_dotw(i) = (-cos(psi)*cos(theta)/sin(theta))*w(1) + (-sin(psi)*cos(theta)/sin(theta))*w(2) + w(3);

                M(:,1) = R0*m(:,i)./norm(m(:,i));
                M(:,2) = R1*m(:,i+1)./norm(m(:,i+1));

                psi1 = atan2(M(2,1),M(1,1));
                psi2 = atan2(M(2,2),M(1,2));
                psi_dotm(i) = (psi2-psi1)/dt;
            end
        end

        if sum(abs(psi_dotw - psi_dotm) > 10*(pi/180)) ~= 0
            res = 1;
        else
            res = 0;
        end
    end
else
    res = 1;
end


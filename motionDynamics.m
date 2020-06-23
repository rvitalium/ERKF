function res = motionDynamics(x,u,g,s,xind,uind)

%{

possible values for res

0 -> static
1 -> not static: acceleratometer, not gyro
2 -> not static: not accelerometer, gyro
3 -> not static: accelerometer + gyro

%}

a = u(uind.a,:);
w = u(uind.w,:);
wb = x(xind.wb,:);

alpha_ = 0.05;
res = [];

% check acceleration criteria

epsilon_a = 0.2;
for i = 1:size(a,2)
    a_norm(i,1) = sqrt(a(:,i)'*a(:,i));
end
g_norm = norm(g);

if sum(abs(a_norm-g_norm) > epsilon_a) ~= 0
    res = 1;
end

% check angular velocity criteria

% epsilon_w = 0.5;
% for i = size(w,2)
%     w_norm(i,1) = sqrt((w(:,i))'*(w(:,i)));
% end

maDist2 = (w-wb)'/(s^2.*eye(3))*(w-wb);
r = maDist2 < chi2inv(1-alpha_,size(maDist2,2));

if r == 0
    if isempty(res)
        res = 2;
    else
        res = res + 2;
    end
end

if isempty(res)
    res = 0;
end

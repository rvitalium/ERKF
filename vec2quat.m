function q = vec2quat(omega,dt)
    phi = norm(omega)*dt;
    u = omega./norm(omega);
    q = [cos(phi/2); u.*sin(phi/2)];
end
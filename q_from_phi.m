function [quat] = q_from_phi(phi)

phi = phi(:);

theta = norm(phi);
r = phi/theta;

quat = [r*sin(theta/2);cos(theta/2)];

if theta == 0 
    quat = [0 0 0 1]';
end

end

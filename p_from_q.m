function [p] = p_from_q(q)
% convert a quaternion, scalar last, to modified rodgrigues parameter

q = q(:);
v = q(1:3);
s = q(4);

p = v/(1+s);

if norm(p)>1
    p = -p/(norm(p)^2);
end
end

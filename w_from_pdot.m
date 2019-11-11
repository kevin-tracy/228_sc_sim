function [w] = w_from_pdot(p,pdot)
% this is the kinematics of the modified rodrigues parameter assuming that
% attitude is being denoted as N_R_B using the kane/levinson convention
p = p(:);
pdot = pdot(:);

w = (4/(1+norm(p)^2))*(eye(3) + 2*(hat(p)^2 - hat(p))/(1+norm(p)^2))*pdot;

end
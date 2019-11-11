function DCM = dcm_from_q(quat)
%this function takes a unit quaternion SCALAR LAST and gives the DCM


quat = quat(:);

v = quat(1:3);
s = quat(4);

DCM = eye(3) +2*hat(v)*(s*eye(3) + hat(v));

end

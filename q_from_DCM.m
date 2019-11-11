function [q] = q_from_DCM(R)
% Scalar Last, Kane/Levinson Convention
T = R(1,1) + R(2,2) + R(3,3);
if T> R(1,1) && T> R(2,2) && T>R(3,3)
    q4 = .5*sqrt(1+T);
    r  = .25/q4;
    q1 = (R(3,2) - R(2,3))*r;
    q2 = (R(1,3) - R(3,1))*r;
    q3 = (R(2,1) - R(1,2))*r;
elseif R(1,1)>R(2,2) && R(1,1)>R(3,3)
    q1 = .5*sqrt(1-T + 2*R(1,1));
    r  = .25/q1;
    q4 = (R(3,2) - R(2,3))*r;
    q2 = (R(1,2) + R(2,1))*r;
    q3 = (R(1,3) + R(3,1))*r;
elseif R(2,2)>R(3,3)
    q2 = .5*sqrt(1-T + 2*R(2,2));
    r  = .25/q2;
    q4 = (R(1,3) - R(3,1))*r;
    q1 = (R(1,2) + R(2,1))*r;
    q3 = (R(2,3) + R(3,2))*r;
else
    q3 = .5*sqrt(1-T + 2*R(3,3));
    r  = .25/q3;
    q4 = (R(2,1) - R(1,2))*r;
    q1 = (R(1,3) + R(3,1))*r;
    q2 = (R(2,3) + R(3,2))*r;
end
q = [q1;q2;q3;q4];
if q4<0
    q = -q;
end
end
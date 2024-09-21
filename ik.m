function [q_fb] = ik(lm,qn,Ln)

% Inverse kinematics of full body

lm = [-lm(:,3), lm(:,1), -lm(:,2)]; % Fixing of the coordinates
lm = lm(2:end,:);

%% Calculation of the first two degrees of freedom and rotation matrices

ltd = (lm(23,:) + lm(24,:))/2;
lm = lm - ltd;

yB = (lm(23,:) - lm(24,:))/norm(lm(23,:)-lm(24,:));
qb = [atan(-yB(1)/yB(2)); asin(yB(3))];

R1 = [cos(qb(1)) -sin(qb(1)) 0;
    sin(qb(1)) cos(qb(1)) 0;
    0 0 1];
R2 = [1 0 0;
    0 cos(qb(2)) -sin(qb(2));
    0 sin(qb(2)) cos(qb(2))];
R = R1*R2;
lm = lm*R;

%% Upper body

xgn = [lm(12,:), lm(14,:), lm(16,:), lm(11,:), lm(13,:), lm(15,:)].';
lbound_w = -pi/6;
ubound_w = pi/4;
lbound_a = [-pi/2;-3*pi/4;0];
ubound_a = [pi;pi/4;5*pi/6];
lbound_ub = [-Inf; lbound_w; lbound_a; -ubound_a];
ubound_ub = [Inf; ubound_w; ubound_a; -lbound_a];

q_ub = ik_ub(Ln(1:5),qn(3:10,1),xgn,lbound_ub,ubound_ub);

%% Rotation according to third degree of freedom, to align with hips

R3 = [cos(q_ub(1)) 0 sin(q_ub(1));
    0 1 0;
    -sin(q_ub(1)) 0 cos(q_ub(1))];
lm = lm*R3;

%% Right leg

xgn = [lm(26,:), lm(28,:), lm(32,:)].';
lbound_l = [-pi/4;-pi/3;-3*pi/4;-pi/4];
ubound_l = [3*pi/4;pi/6;0;pi/6];
q_rl = ik_rl(Ln(7:10),qn(11:14,1),xgn,lbound_l,ubound_l);

%% Left leg

xgn = [lm(25,:), lm(27,:), lm(31,:)].';
q_ll = ik_ll(Ln(7:10),qn(15:18),xgn,-ubound_l,-lbound_l);

%% 
q_fb = [qb; q_ub; q_rl; q_ll];
toc

end
function [q_fb] = ik(lm,Ln)

lm = [-lm(:,3), lm(:,1), -lm(:,2)]; % popravka koordinata

lm = lm(2:end,:);

%%

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

%%

xgn = [lm(12,:), lm(14,:), lm(16,:), (lm(18,:)+lm(20,:))/2, lm(11,:), lm(13,:), lm(15,:), (lm(17,:)+lm(19,:))/2].';
qn = zeros(18,1);
lbound_w = [-pi/3;-pi/6;-pi/2];
ubound_w = [pi/3;pi/2;pi/2];
lbound_a = [-pi/2;-3*pi/4;-pi/2;0;-pi/2;-pi/2;-pi/6];
ubound_a = [pi;pi/4;pi/2;5*pi/6;pi/2;pi/2;pi/6];
lbound_ub = [-Inf; lbound_w; lbound_a; -ubound_a];
ubound_ub = [Inf; ubound_w; ubound_a; -lbound_a];

q_ub = ik_ub(Ln(1:6),qn,xgn,lbound_ub,ubound_ub);

%%
R3 = [cos(q_ub(1)) 0 sin(q_ub(1));
    0 1 0;
    -sin(q_ub(1)) 0 cos(q_ub(1))];
lm = lm*R3;

%%

xgn = [lm(26,:), lm(28,:), lm(30,:), lm(32,:)].';
qn = zeros(7,1);
lbound_l = [-pi/4;-pi/3;-pi/2;-3*pi/4;-pi/4;-pi/6;-pi/6];
ubound_l = [3*pi/4;pi/6;pi/2;0;pi/6;pi/6;pi/6];

q_rl = ik_rl(Ln(7:13),qn,xgn,lbound_l,ubound_l);


%%
xgn = [lm(25,:), lm(27,:), lm(29,:), lm(31,:)].';
qn = zeros(7,1);

q_ll = ik_ll(Ln(7:13),qn,xgn,-ubound_l,-lbound_l);


%%
q_fb = [qb; q_ub; q_rl; q_ll];
toc
end
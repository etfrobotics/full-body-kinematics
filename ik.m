function [q_fb] = ik(lm,qn,Ln)

lm = [-lm(:,3), lm(:,1), -lm(:,2)]; % popravka koordinata

lm = lm(2:end,:);

% Lbt = norm((lm(11,:)+lm(12,:))/2 - (lm(23,:)+lm(24,:))/2);
% Lbw = Lbt*0.5; Lwt = Lbt*0.5;
% Lta = norm(lm(11,:)-lm(12,:))/2;
% La1 = (norm(lm(12,:)-lm(14,:)) + norm(lm(11,:)-lm(13,:)))/2;
% La2 = (norm(lm(14,:)-lm(16,:)) + norm(lm(13,:)-lm(15,:)))/2;
% 
% Lbl = norm(lm(23,:)-lm(24,:))/2;
% Ll1 = (norm(lm(23,:)-lm(25,:)) + norm(lm(24,:)-lm(26,:)))/2;
% Ll2 = (norm(lm(25,:)-lm(29,:)) + norm(lm(26,:)-lm(30,:)))/2;
% Ll3 = (norm(lm(29,:)-lm(31,:)) + norm(lm(30,:)-lm(32,:)))/2;
% 
% Ln = [Lbw;Lwt;Lta;La1;La2;Lbl;Ll1;Ll2;Ll3];


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

xgn = [lm(12,:), lm(14,:), lm(16,:), lm(11,:), lm(13,:), lm(15,:)].';
lbound_w = -pi/6;
ubound_w = pi/2;
lbound_a = [-pi/2;-3*pi/4;0];
ubound_a = [pi;pi/4;5*pi/6];
lbound_ub = [-Inf; lbound_w; lbound_a; -ubound_a];
ubound_ub = [Inf; ubound_w; ubound_a; -lbound_a];

q_ub = ik_ub(Ln(1:5),qn(3:10,1),xgn,lbound_ub,ubound_ub);

%%
R3 = [cos(q_ub(1)) 0 sin(q_ub(1));
    0 1 0;
    -sin(q_ub(1)) 0 cos(q_ub(1))];
lm = lm*R3;

%%

xgn = [lm(26,:), lm(30,:), lm(32,:)].';
lbound_l = [-pi/4;-pi/3;-3*pi/4;-pi/4];
ubound_l = [3*pi/4;pi/6;0;pi/6];

% Ll1 = norm(lm(24,:)-lm(26,:));
% Ll2 = norm(lm(26,:)-lm(30,:));
% Ll3 = norm(lm(30,:)-lm(32,:));
% Ln = [Lbw;Lwt;Lta;La1;La2;Lbl;Ll1;Ll2;Ll3];

q_rl = ik_rl(Ln(7:10),qn(11:14,1),xgn,lbound_l,ubound_l);


%%
xgn = [lm(25,:), lm(29,:), lm(31,:)].';

% Ll1 = norm(lm(23,:)-lm(25,:));
% Ll2 = norm(lm(25,:)-lm(29,:));
% Ll3 = norm(lm(29,:)-lm(31,:));
% Ln = [Lbw;Lwt;Lta;La1;La2;Lbl;Ll1;Ll2;Ll3];

q_ll = ik_ll(Ln(7:10),qn(15:18),xgn,-ubound_l,-lbound_l);


%%
q_fb = [qb; q_ub; q_rl; q_ll];
toc
end
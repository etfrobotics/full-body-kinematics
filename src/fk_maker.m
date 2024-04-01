clc
clear
close all

%% Base frame

q = sym('q',[18,1]);
L = sym('L',[6,1]);
R3 = [cos(q(1)) 0 sin(q(1));
    0 1 0;
    -sin(q(1)) 0 cos(q(1))];

H_b = sym(zeros(4,4));
H_b(1:3,1:3) = R3;
H_b(4,4) = 1;

%% Waist 3DOF

H_w = sym(zeros(4,4,3));
H_w(:,:,1) = H_b*dh_calc(0,L(1),0,pi/2,pi/2);
H_w(:,:,2) = H_w(:,:,1)*dh_calc(q(2),0,0,pi/2,pi/2);
H_w(:,:,3) = H_w(:,:,2)*dh_calc(q(3),0,0,pi/2,pi/2);
H_t = H_w(:,:,3)*dh_calc(q(4),L(2),0,0,0);


%% Right arm 7DOF

d_ra = [0 0 L(4) 0 L(5) 0 0];
a_ra = [0 0 0 0 0 0 L(6)];
alpha_ra = [pi/2 -pi/2 -pi/2 pi/2 pi/2 pi/2 0];
offset_ra = [pi/2 pi/2 -pi/2 0 pi/2 pi/2 0];

Htemp = [1,0,0,0;
        0,0,-1,-L(3);
        0,1,0,0;
        0,0,0,1];
H_ra = sym(zeros(4,4,8));
H_ra(:,:,1) = H_t*Htemp;
for i = 1:7
    H_ra(:,:,i+1) = H_ra(:,:,i)*dh_calc(q(i+4), d_ra(i), a_ra(i), alpha_ra(i), offset_ra(i));
end

%% Left arm 7DOF

d_la = [0 0 L(4) 0 L(5) 0 0];
a_la = [0 0 0 0 0 0 L(6)];
alpha_la = -alpha_ra;
offset_la = -offset_ra;

Htemp = [1,0,0,0;
        0,0,1,L(3);
        0,-1,0,0;
        0,0,0,1];
H_la = sym(zeros(4,4,8));
H_la(:,:,1) = H_t*Htemp;
for i = 1:7
    H_la(:,:,i+1) = H_la(:,:,i)*dh_calc(q(i+11), d_la(i), a_la(i), alpha_la(i), offset_la(i));
end

%%
tic
fk_ub = cat(3,H_b,H_w,H_t,H_ra,H_la);
toc
tic
matlabFunction(fk_ub, 'File', 'fk_ub_computable.m', 'Vars', {q, L});
toc
%%
xg = sym('xg', [24, 1]);
J_ub = sum(([fk_ub(1:3,4,6); fk_ub(1:3,4,9); fk_ub(1:3,4,11); fk_ub(1:3,4,13); ...
    fk_ub(1:3,4,14); fk_ub(1:3,4,17); fk_ub(1:3,4,19); fk_ub(1:3,4,21)] - xg).^2);
tic
dJ_ub_dq = jacobian(J_ub, q);
toc
tic
ddJ_ub_ddq = hessian(J_ub, q);
toc
%%
tic
matlabFunction(J_ub, 'File', 'J_ub_computable.m', 'Vars', {q, L, xg});
toc
tic
matlabFunction(dJ_ub_dq, 'File', 'dJ_ub_dq_computable.m', 'Vars', {q, L, xg});
toc
% tic
% matlabFunction(ddJ_ub_ddq, 'File', 'ddJ_ub_ddq_computable.m', 'Vars', {q, L, xg});
% toc

%% Right leg 7DOF

q_l = sym('ql', [7, 1]);
L_l= sym('Ll', [7, 1]);

d_rl = [0 0 L_l(2) 0 0 -L_l(5)];
a_rl = [0 0 0 L_l(3) 0 L_l(4)];
alpha_rl = [pi/2 -pi/2 -pi/2 0 -pi/2 0];
offset_rl = [pi/2 pi/2 -pi/2 -pi/2 0 0];

Htemp = [1,0,0,0;
        0,0,-1,-L_l(1);
        0,1,0,0;
        0,0,0,1];
H_rl = sym(zeros(4,4,9));
H_rl(:,:,1) = Htemp;
for i = 1:6
    H_rl(:,:,i+1) = H_rl(:,:,i)*dh_calc(q_l(i), d_rl(i), a_rl(i), alpha_rl(i), offset_rl(i));
end

H_rl(:,:,8) = H_rl(:,:,7)*dh_calc(0,L_l(5)+L_l(6),0,pi/2,0);
H_rl(:,:,9) = H_rl(:,:,8)*dh_calc(q_l(7),0,L_l(7),-pi/2,pi/2);

fk_rl = H_rl;

tic
matlabFunction(fk_rl, 'File', 'fk_rl_computable.m', 'Vars', {q_l, L_l});
toc
%%
xg = sym('xg', [12, 1]);
J_rl = sum(([fk_rl(1:3,4,4); fk_rl(1:3,4,5); fk_rl(1:3,4,7); fk_rl(1:3,4,9)] - xg).^2);
tic
dJ_rl_dq = jacobian(J_rl, q_l);
toc
tic
ddJ_rl_ddq = hessian(J_rl, q_l);
toc
%%
tic
matlabFunction(J_rl, 'File', 'J_rl_computable.m', 'Vars', {q_l, L_l, xg});
toc
tic
matlabFunction(dJ_rl_dq, 'File', 'dJ_rl_dq_computable.m', 'Vars', {q_l, L_l, xg});
toc
tic
matlabFunction(ddJ_rl_ddq, 'File', 'ddJ_rl_ddq_computable.m', 'Vars', {q_l, L_l, xg});
toc

%% Left leg 7DOF

q_l = sym('ql', [7, 1]);
L_l= sym('Ll', [7, 1]);

d_ll = [0 0 L_l(2) 0 0 -L_l(5)];
a_ll = [0 0 0 L_l(3) 0 L_l(4)];
alpha_ll = -alpha_rl;
offset_ll = -offset_rl;

Htemp = [1,0,0,0;
        0,0,1,L_l(1);
        0,-1,0,0;
        0,0,0,1];
H_ll = sym(zeros(4,4,9));
H_ll(:,:,1) = Htemp;
for i = 1:6
    H_ll(:,:,i+1) = H_ll(:,:,i)*dh_calc(q_l(i), d_ll(i), a_ll(i), alpha_ll(i), offset_ll(i));
end

H_ll(:,:,8) = H_ll(:,:,7)*dh_calc(0,L_l(5)+L_l(6),0,-pi/2,0);
H_ll(:,:,9) = H_ll(:,:,8)*dh_calc(q_l(7),0,L_l(7),pi/2,-pi/2);

fk_ll = H_ll;

tic
matlabFunction(fk_ll, 'File', 'fk_ll_computable.m', 'Vars', {q_l, L_l});
toc
%%
xg = sym('xg', [12, 1]);
J_ll = sum(([fk_ll(1:3,4,4); fk_ll(1:3,4,5); fk_ll(1:3,4,7); fk_ll(1:3,4,9)] - xg).^2);
tic
dJ_ll_dq = jacobian(J_ll, q_l);
toc
tic
ddJ_ll_ddq = hessian(J_ll, q_l);
toc
%%
tic
matlabFunction(J_ll, 'File', 'J_ll_computable.m', 'Vars', {q_l, L_l, xg});
toc
tic
matlabFunction(dJ_ll_dq, 'File', 'dJ_ll_dq_computable.m', 'Vars', {q_l, L_l, xg});
toc
tic
matlabFunction(ddJ_ll_ddq, 'File', 'ddJ_ll_ddq_computable.m', 'Vars', {q_l, L_l, xg});
toc

function H = fk(q,L)
sz = size(q);
if sz(2) > sz(1)
    q = q.';
end

R1 = [cos(q(1)) -sin(q(1)) 0;
    sin(q(1)) cos(q(1)) 0;
    0 0 1];
R2 = [1 0 0;
    0 cos(q(2)) -sin(q(2));
    0 sin(q(2)) cos(q(2))];
R = R1*R2;
fk_ub = fk_ub_computable(q(3:20),L(1:6));

Htemp = eye(4);
Htemp(1:3,1:3) = R;
for i = 1:length(fk_ub)
    fk_ub(:,:,i) = Htemp*fk_ub(:,:,i);
end

R3 = [cos(q(3)) 0 sin(q(3));
    0 1 0;
    -sin(q(3)) 0 cos(q(3))];
R = R1*R2*R3;

fk_rl = fk_rl_computable(q(21:27),L(7:13));

Htemp = eye(4);
Htemp(1:3,1:3) = R;
for i = 1:length(fk_rl)
    fk_rl(:,:,i) = Htemp*fk_rl(:,:,i);
end

fk_ll = fk_ll_computable(q(28:34),L(7:13));

Htemp = eye(4);
Htemp(1:3,1:3) = R;
for i = 1:length(fk_ll)
    fk_ll(:,:,i) = Htemp*fk_ll(:,:,i);
end

H = cat(3,fk_ub,fk_rl,fk_ll);

end
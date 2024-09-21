function H_CoM = fk_CoM(H,q,L,BSIP,M)

% Function for calculating CoM positions based on joint positions

d = BSIP.torso_head.r_com;
H_temp = eye(4); H_temp(1:3,4) = d;
CoM_torso_head = H(:,:,3)*H_temp;
Rot = [1 0 0; 0 0 -1; 0 1 0]; d = BSIP.pelvis_abdomen.r_com;
H_temp = eye(4); H_temp(1:3,1:3) = Rot; H_temp(1:3,4) = Rot*d;
CoM_pelvis_abdomen = H(:,:,2)*H_temp;

Rot = [0 0 -1; 1 0 0; 0 -1 0]; d = BSIP.r_upper_arm.r_com+[0;0;L(4)];
H_temp = eye(4); H_temp(1:3,1:3) = Rot; H_temp(1:3,4) = Rot*d;
CoM_r_upper_arm = H(:,:,6)*H_temp;
Rot = [0 0 -1; 1 0 0; 0 -1 0]; d = BSIP.r_forearm.r_com+[0;0;L(5)];
H_temp = eye(4); H_temp(1:3,1:3) = Rot; H_temp(1:3,4) = Rot*d;
CoM_r_forearm = H(:,:,7)*H_temp;

Rot = [0 0 -1; -1 0 0; 0 1 0]; d = BSIP.l_upper_arm.r_com+[0;0;L(4)];
H_temp = eye(4); H_temp(1:3,1:3) = Rot; H_temp(1:3,4) = Rot*d;
CoM_l_upper_arm = H(:,:,10)*H_temp;
Rot = [0 0 -1; -1 0 0; 0 1 0]; d = BSIP.l_forearm.r_com+[0;0;L(5)];
H_temp = eye(4); H_temp(1:3,1:3) = Rot; H_temp(1:3,4) = Rot*d;
CoM_l_forearm = H(:,:,11)*H_temp;

Rot = [0 0 -1; 1 0 0; 0 -1 0]; d = BSIP.r_thigh.r_com+[0;0;L(8)];
H_temp = eye(4); H_temp(1:3,1:3) = Rot; H_temp(1:3,4) = Rot*d;
CoM_r_thigh = H(:,:,14)*H_temp;
Rot = [0 0 -1; 1 0 0; 0 -1 0]; d = BSIP.r_shank.r_com+[0;0;L(9)];
H_temp = eye(4); H_temp(1:3,1:3) = Rot; H_temp(1:3,4) = Rot*d;
CoM_r_shank = H(:,:,15)*H_temp;
Rot = [0 0 -1; 1 0 0; 0 -1 0]*[cos(-q(14)) 0 sin(-q(14)); 0 1 0; -sin(-q(14)) 0 cos(-q(14))]; d = BSIP.r_foot.r_com;
H_temp = eye(4); H_temp(1:3,1:3) = Rot; H_temp(1:3,4) = Rot*d;
CoM_r_foot = H(:,:,15)*H_temp;

Rot = [0 0 -1; -1 0 0; 0 1 0]; d = BSIP.l_thigh.r_com+[0;0;L(8)];
H_temp = eye(4); H_temp(1:3,1:3) = Rot; H_temp(1:3,4) = Rot*d;
CoM_l_thigh = H(:,:,19)*H_temp;
Rot = [0 0 -1; -1 0 0; 0 1 0]; d = BSIP.l_shank.r_com+[0;0;L(9)];
H_temp = eye(4); H_temp(1:3,1:3) = Rot; H_temp(1:3,4) = Rot*d;
CoM_l_shank = H(:,:,20)*H_temp;
Rot = [0 0 -1; -1 0 0; 0 1 0]*[cos(q(18)) 0 sin(q(18)); 0 1 0; -sin(q(18)) 0 cos(q(18))]; d = BSIP.l_foot.r_com;
H_temp = eye(4); H_temp(1:3,1:3) = Rot; H_temp(1:3,4) = Rot*d;
CoM_l_foot = H(:,:,20)*H_temp;


H_CoM = cat(3,CoM_pelvis_abdomen,CoM_torso_head, ...
    CoM_r_upper_arm,CoM_r_forearm,CoM_l_upper_arm,CoM_l_forearm, ...
    CoM_r_thigh,CoM_r_shank,CoM_r_foot,CoM_l_thigh,CoM_l_shank,CoM_l_foot);

fn = fieldnames(BSIP);
d_CoM = [0;0;0];

for i = 1:length(fn)
    d_CoM = d_CoM + BSIP.(fn{i}).m*H_CoM(1:3,4,i)/M;
end
H_temp = zeros(4,4); H_temp(1:3,4) = d_CoM;
CoM = H(:,:,1) + H_temp;

H_CoM = cat(3,H_CoM,CoM);
   
end
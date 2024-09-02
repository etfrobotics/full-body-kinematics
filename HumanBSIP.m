function bsipStruct = HumanBSIP(M,Ln)

% Function for calculating body segment inertial parameters based on 
% "Estimation of the Body Segment Inertial Parameters for the Rigid Body 
% Biomechanical Models Used in Motion Analysis" by R. Dumas and J. Wojtusch
% 
% The function calculates body segment weights, relative positions of the
% centers of mass and matrices of inertia

segmentNames = [    % Ordered as in Dumas (2018)
    "head";
    "torso_head";
    "abdomen";
    "pelvis_abdomen";
    "r_upper_arm";
    "r_forearm";
    "r_hand";
    "r_thigh";
    "r_shank";
    "r_foot";
    "l_upper_arm";
    "l_forearm";
    "l_hand";
    "l_thigh";
    "l_shank";
    "l_foot";
    ];

% Initialize structure in correct order
bsipStruct.pelvis_abdomen.m = 0; bsipStruct.torso_head.m = 0; 
bsipStruct.r_upper_arm.m = 0; bsipStruct.r_forearm.m = 0;
bsipStruct.l_upper_arm.m = 0; bsipStruct.l_forearm.m = 0;
bsipStruct.r_thigh.m = 0; bsipStruct.r_shank.m = 0; bsipStruct.r_foot.m = 0;
bsipStruct.l_thigh.m = 0; bsipStruct.l_shank.m = 0; bsipStruct.l_foot.m = 0;

L_head=0.8*Ln(2); L_torso=Ln(2); L_abdomen=0.6*Ln(1); L_pelvis=0.4*Ln(1); L_upper_arm=Ln(4);
L_forearm=Ln(5); L_hand=Ln(6); L_thigh=Ln(8); L_shank=Ln(9); L_foot=Ln(10);

L = [L_head; L_torso; L_abdomen; L_pelvis; L_upper_arm; L_forearm; L_hand; L_thigh; L_shank; L_foot];
p = [6.7; 30.4; 2.9; 14.2; 2.4; 1.7; 0.6; 12.3; 4.8; 1.2]/100;
cx = [2; 0; 17.6; -0.2; 1.8; -1.3; 8.2; -4.1; -4.8; 50.2]/100;
cy = [53.4; -55.5; -36.1; -28.2; -48.2; -41.7; -83.9; -42.9; -41; -19.9]/100;
cz = [0.1; -0.4; -3.3; -0.6; -3.1; 1.1; 7.5; 3.3; 0.7; 3.4]/100;
rxx = [28; 42; 54; 102; 29; 28; 61; 29; 28; 22]/100;
ryy = [21; 33; 66; 106; 13; 11; 38; 15; 10; 49]/100;
rzz = [30; 36; 40; 96; 30; 28; 56; 30; 28; 48]/100;
rxy = [7*1i; 11*1i; 11; 25*1i; 5; 8; 22; 7; 4*1i; 17]/100;
rxz = [2*1i; 1; 6*1i; 12*1i; 3; 1*1i; 15; 2*1i; 2*1i; 11*1i]/100;
ryz = [3; 3; 5*1i; 8*1i; 13*1i; 2; 20*1i; 7*1i; 4; 0]/100;

P = [1 0 0; 0 0 -1; 0 1 0];
for i = 1:10
    bsipStruct.(segmentNames(i)).m = p(i)*M;
    bsipStruct.(segmentNames(i)).r_com = P*L(i)*[cx(i); cy(i); cz(i)];
    bsipStruct.(segmentNames(i)).I = bsipStruct.(segmentNames(i)).m*L(i).^2*...
        P*[rxx(i)^2, rxy(i)^2, rxz(i)^2;
        rxy(i)^2, ryy(i)^2, ryz(i)^2;
        rxz(i)^2, ryz(i)^2, rzz(i)^2;]*P.';
end

S = [1 0 0; 0 -1 0; 0 0 1];
for i = 11:16
    bsipStruct.(segmentNames(i)).m = p(i-6)*M;
    bsipStruct.(segmentNames(i)).r_com = S*bsipStruct.(segmentNames(i-6)).r_com;
    bsipStruct.(segmentNames(i)).I = S*bsipStruct.(segmentNames(i-6)).I*S.';  
end

% Joining of certain body segments due to omitted degrees of freedom
torso_head_com = (bsipStruct.("head").r_com.*bsipStruct.("head").m + ...
    bsipStruct.("torso_head").r_com.*bsipStruct.("torso_head").m)/(bsipStruct.("head").m+bsipStruct.("torso_head").m);
r_th_head = bsipStruct.("head").r_com-torso_head_com; r_th_torso = bsipStruct.("torso_head").r_com-torso_head_com;
bsipStruct.("torso_head").I = bsipStruct.("head").I + bsipStruct.("torso_head").I - ...
    bsipStruct.("head").m.*cross([r_th_head r_th_head r_th_head],eye(3))^2 - ...
    bsipStruct.("torso_head").m.*cross([r_th_torso r_th_torso r_th_torso],eye(3))^2; 
bsipStruct.("torso_head").m = bsipStruct.("head").m + bsipStruct.("torso_head").m;
bsipStruct.("torso_head").r_com = torso_head_com;

bsipStruct.("pelvis_abdomen").r_com = bsipStruct.("pelvis_abdomen").r_com -[0;0;L_abdomen];
pelvis_abdomen_com = (bsipStruct.("abdomen").r_com.*bsipStruct.("abdomen").m + ...
    bsipStruct.("pelvis_abdomen").r_com.*bsipStruct.("pelvis_abdomen").m)/(bsipStruct.("abdomen").m+bsipStruct.("pelvis_abdomen").m);
r_pa_abdomen = bsipStruct.("abdomen").r_com-pelvis_abdomen_com; r_pa_pelvis = bsipStruct.("pelvis_abdomen").r_com-pelvis_abdomen_com;
bsipStruct.("pelvis_abdomen").I = bsipStruct.("abdomen").I + bsipStruct.("pelvis_abdomen").I - ...
    bsipStruct.("abdomen").m.*cross([r_pa_abdomen r_pa_abdomen r_pa_abdomen],eye(3))^2 - ...
    bsipStruct.("pelvis_abdomen").m.*cross([r_pa_pelvis r_pa_pelvis r_pa_pelvis],eye(3))^2; 
bsipStruct.("pelvis_abdomen").m = bsipStruct.("abdomen").m + bsipStruct.("pelvis_abdomen").m;
bsipStruct.("pelvis_abdomen").r_com = pelvis_abdomen_com;

bsipStruct.("r_hand").r_com = bsipStruct.("r_hand").r_com -[0;0;L_forearm];
r_forearm_com = (bsipStruct.("r_hand").r_com.*bsipStruct.("r_hand").m + ...
    bsipStruct.("r_forearm").r_com.*bsipStruct.("r_forearm").m)/(bsipStruct.("r_hand").m+bsipStruct.("r_forearm").m);
r_fh_hand = bsipStruct.("r_hand").r_com-r_forearm_com; r_fh_forearm = bsipStruct.("r_forearm").r_com-r_forearm_com;
bsipStruct.("r_forearm").I = bsipStruct.("r_hand").I + bsipStruct.("r_forearm").I - ...
    bsipStruct.("r_hand").m.*cross([r_fh_hand r_fh_hand r_fh_hand],eye(3))^2 - ...
    bsipStruct.("r_forearm").m.*cross([r_fh_forearm r_fh_forearm r_fh_forearm],eye(3))^2; 
bsipStruct.("r_forearm").m = bsipStruct.("r_hand").m + bsipStruct.("r_forearm").m;
bsipStruct.("r_forearm").r_com = r_forearm_com;

bsipStruct.("l_hand").r_com = bsipStruct.("l_hand").r_com -[0;0;L_forearm];
r_forearm_com = (bsipStruct.("l_hand").r_com.*bsipStruct.("l_hand").m + ...
    bsipStruct.("l_forearm").r_com.*bsipStruct.("l_forearm").m)/(bsipStruct.("l_hand").m+bsipStruct.("l_forearm").m);
r_fh_hand = bsipStruct.("l_hand").r_com-r_forearm_com; r_fh_forearm = bsipStruct.("l_forearm").r_com-r_forearm_com;
bsipStruct.("l_forearm").I = bsipStruct.("l_hand").I + bsipStruct.("l_forearm").I - ...
    bsipStruct.("l_hand").m.*cross([r_fh_hand r_fh_hand r_fh_hand],eye(3))^2 - ...
    bsipStruct.("l_forearm").m.*cross([r_fh_forearm r_fh_forearm r_fh_forearm],eye(3))^2; 
bsipStruct.("l_forearm").m = bsipStruct.("l_hand").m + bsipStruct.("l_forearm").m;
bsipStruct.("l_forearm").r_com = r_forearm_com;

bsipStruct = rmfield(bsipStruct,"head");
bsipStruct = rmfield(bsipStruct,"abdomen");
bsipStruct = rmfield(bsipStruct,"r_hand");
bsipStruct = rmfield(bsipStruct,"l_hand");
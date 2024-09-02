clear
close all
clc
    
%% Smoothing of the landmarks

addpath('src')
tic
lm = readmatrix('landmarks.csv');  % Loading of the data

for i = 1:33
    for j = 3:5
        lm(i:33:33*lm(end,1)+i,j) = smoothdata(lm(i:33:33*lm(end,1)+i,j),"gaussian",18);
    end
end

%% Inverse kinematics

% Segment lengths are calculated based on the first frame of the video
lmT = [-lm(2:33,5), lm(2:33,3), -lm(2:33,4)];
Lbt = norm((lmT(11,:)+lmT(12,:))/2 - (lmT(23,:)+lmT(24,:))/2);
Lbw = Lbt*0.5; Lwt = Lbt*0.5;
Lta = norm(lmT(11,:)-lmT(12,:))/2;
La1 = (norm(lmT(12,:)-lmT(14,:)) + norm(lmT(11,:)-lmT(13,:)))/2;
La2 = (norm(lmT(14,:)-lmT(16,:)) + norm(lmT(13,:)-lmT(15,:)))/2;
Lh = (norm(lmT(16,:)-lmT(20,:)) + norm(lmT(15,:)-lmT(19,:)))/2;

Lbl = norm(lmT(23,:)-lmT(24,:))/2;
Ll1 = (norm(lmT(23,:)-lmT(25,:)) + norm(lmT(24,:)-lmT(26,:)))/2;
Ll2 = (norm(lmT(25,:)-lmT(29,:)) + norm(lmT(26,:)-lmT(30,:)))/2;
Ll3 = (norm(lmT(29,:)-lmT(31,:)) + norm(lmT(30,:)-lmT(32,:)))/2;

Ln = [Lbw;Lwt;Lta;La1;La2;Lh;Lbl;Ll1;Ll2;Ll3];
qn = zeros(18,1);
bsip = HumanBSIP(80,Ln);
q = zeros(lm(end,1)+1,18);
r_CoM = zeros(lm(end,1)+1,3);

for i = 0:lm(end,1)
    q(i+1,:) = ik(lm(i*33+1:(i+1)*33,3:5),qn,Ln);
    qn = q(i+1,:).';
    H_CoM = fk_CoM(fk(qn,Ln),Ln,bsip,80);
    n = length(H_CoM);
    r_CoM(i+1,:) = H_CoM(1:3,4,13)-[-lm(i*33+32,5)-lm(i*33+33,5); lm(i*33+32,3)+lm(i*33+33,3); -lm(i*33+32,4)-lm(i*33+33,4)]/2;
end

% plotAngles(q)

%% Filtering of the angles

q = cat(1,repmat(q(1,:),50,1),q,repmat(q(end,:),50,1));
r_CoM = cat(1,repmat(r_CoM(1,:),50,1),r_CoM,repmat(r_CoM(end,:),50,1));
for i = 1:18
    q_s(:,i) = smoothdata(q(:,i),"gaussian",9);
end

for i = 1:18
    q_lp(:,i) = lowpass(q(:,i),2,90);
end
for i = 1:3
    r_CoM(:,i) = lowpass(r_CoM(:,i),2,90);
end

q = q_lp(51:end-50,:);
r_CoM = r_CoM(51:end-50,:);

%% Comparison of left and right angles

plotAngles(q)
toc

%knee f-e
figure()
plot(-q(:,13))
hold on
plot(q(:,17))

%hip f-e
figure()
plot(q(:,11))
hold on
plot(-q(:,15))

%hip a-a
figure()
plot(-q(:,12))
hold on
plot(q(:,16))

%% Differentiaton of angles and CoM

t = linspace(0,1/90*size(q,1),size(q,1)).';
for i=1:18
 qd(:,i) = diff(q(:,i)) ./diff(t);
 qdd(:,i) = diff(qd(:,i)) ./ diff(t(1:end-1));
end
qd = [zeros(1,18); qd];
qdd = [zeros(2,18); qdd];
for i=1:3
    v_CoM(:,i) = diff(r_CoM(:,i)) ./ diff(t);
end
v_CoM = [zeros(1,3); v_CoM];
v_CoM(:,4:6) = [qd(:,2), qd(:,3), qd(:,1)];
figure()
plot(t(1:end),v_CoM); legend(); hold on

%% Calculating parabolic motion for CoM

vz = [0; v_CoM(:,3)];
z_CoM = r_CoM(:,3);
[vz_max,i_max] = max(vz);
vz(i_max) = 3.1;
for i = i_max+1:181
    vz(i) = vz(i-1)-9.81/90;
    z_CoM(i) = z_CoM(i-1)+vz(i-1)/90;
end

figure()
plot(r_CoM(:,3)); hold on
plot(z_CoM)
plot(vz)

%% Forward kinematics for animation and r_base

close all

vid = [];
videoFileName = fullfile('video.avi');
vw = VideoWriter(videoFileName,'Motion JPEG AVI');
vw.FrameRate = 90;
open(vw);
figure()
displacement = z_CoM - r_CoM(:,3);
for i = 0:lm(end,1)
    qn = q(i+1,:);
    disp = [(lm(i*33+32,5)+lm(i*33+33,5))/2;-(lm(i*33+32,3)+lm(i*33+33,3))/2;(lm(i*33+32,4)+lm(i*33+33,4))/2+displacement(i+1)];
    scatter3(-lm(i*33+2:(i+1)*33,5)+disp(1)*ones(32,1), ...
        lm(i*33+2:(i+1)*33,3)+disp(2)*ones(32,1), ...
        -lm(i*33+2:(i+1)*33,4)+disp(3)*ones(32,1));
    hold on
    H_CoM = fk_CoM(fk(qn,Ln),Ln,bsip,80);
    fk_all = cat(3,fk(q(i+1,:),Ln),H_CoM);
    H_temp = zeros(4,4); H_temp(1:3,4) = [0; 0; displacement(i+1)] + ...
        -[-lm(i*33+32,5)-lm(i*33+33,5); lm(i*33+32,3)+lm(i*33+33,3); -lm(i*33+32,4)-lm(i*33+33,4)]/2;
    CoM = repmat(H_temp,[1 1 length(fk_all)]);
    n = length(H_CoM);
    scatter3(reshape(H_CoM(1,4,1:n),[n,1])+disp(1)*ones(n,1), ...
    reshape(H_CoM(2,4,1:n),[n,1])+disp(2)*ones(n,1), ...
    reshape(H_CoM(3,4,1:n),[n,1])+disp(3)*ones(n,1),'r')
    fk_all = fk_all+CoM;
    plotJoints(fk_all);
    axis([-1 1 -1 1 0 2])
    hold off
    fr = getframe(gcf);
    vid = [vid, fr];
    r_base(i+1,:) = fk_all(1:3,4,1);
    fr = getframe(gcf);
    vid = [vid, fr];
    im = frame2im(fr);

    % Write to the Video File 
    writeVideo(vw,im);
end
close(vw)


%% Calculating v_base in base frame

close all
v_base = [];
for i=1:3
    v_base(:,i) = diff(r_base(:,i)) ./ diff(t);
end
v_base = [zeros(1,3); v_base];

v_base(:,4:6) = v_CoM(:,4:6);

for i = 1:6
    v_base(:,i) = lowpass(v_base(:,i),2,90);
    a_base(:,i) = diff(v_base(:,i)) ./ diff(t);
end
a_base = [zeros(1,6); a_base];
figure()
plot(v_base)
legend
ylim([-5 5])

for i = 0:lm(end,1)
    qn = q(i+1,:);
    disp = [(lm(i*33+32,5)+lm(i*33+33,5))/2;-(lm(i*33+32,3)+lm(i*33+33,3))/2;(lm(i*33+32,4)+lm(i*33+33,4))/2+displacement(i+1)];
    H_CoM = fk_CoM(fk(qn,Ln),Ln,bsip,80);
    fk_all = cat(3,fk(q(i+1,:),Ln),H_CoM);
    H_temp = zeros(4,4); H_temp(1:3,4) = [0; 0; displacement(i+1)] + ...
        -[-lm(i*33+32,5)-lm(i*33+33,5); lm(i*33+32,3)+lm(i*33+33,3); -lm(i*33+32,4)-lm(i*33+33,4)]/2;
    CoM = repmat(H_temp,[1 1 length(fk_all)]);
    fk_all = fk_all+CoM;
    v_base(i+1,1:3) = (fk_all(1:3,1:3,1)*v_base(i+1,1:3).').';
end

figure()
plot(t,v_base(:,1),'r',t,v_base(:,2),'g',t,v_base(:,3),'b',t,v_base(:,4),'r--',t,v_base(:,5),'g--',t,v_base(:,6),'b--','Linewidth',1)
legend('vx','vy','vz','wx','wy','wz')
ylim([-5 5])

%% Defining the body segment graph

lambda = [0 1 2 3 2 5 1 7 8 1 10 11];
joint = [2 4 6 8 10 12 14 15 17 19 20 16 21];
A_graph = zeros(12,13);
for i = 1:12
    for j = 1:11
        if lambda(j+1)==i
            A_graph(i,j) = -1;
        elseif i==j+1
            A_graph(i,j) = 1;
        end
    end
end
A_graph(9,12)=-1; A_graph(12,13)=-1;

figure()
plot(graph(lambda(2:end),2:12))

%% Inverse dynamics

fnames = fieldnames(bsip);
v = zeros(6,12,length(qd));
a = zeros(6,12,length(qd));
f = zeros(72,1);
fint = zeros(78,341);
for i = 1:length(qd) 
    qn = q(i,:);
    qdi = qd(i,:);
    qddi = qdd(i,:);
    % Instead of defining S for each joint
    qdjoint = [0 qdi(4) qdi(5) qdi(7) qdi(8) qdi(10) qdi(11) qdi(13) qdi(14) qdi(15) qdi(17) qdi(18);
        0 0 qdi(6) 0 qdi(9) 0 qdi(12) 0 0 qdi(16) 0 0];
    qddjoint = [0 qddi(4) qddi(5) qddi(7) qddi(8) qddi(10) qddi(11) qddi(13) qddi(14) qddi(15) qddi(17) qddi(18);
        0 0 qddi(6) 0 qddi(9) 0 qddi(12) 0 0 qddi(16) 0 0];
    Sx = [0 0 1 0 1 0 1 0 0 1 0 0];
    Sy = [0 1 -cos(qn(6)) -1 cos(qn(9)) 1 -cos(qn(12)) -1 -1 cos(qn(16)) 1 1];
    Sz = [0 0 sin(qn(6)) 0 -sin(qn(9)) 0 sin(qn(12)) 0 0 -sin(qn(16)) 0 0];
    Stx = [0 0 cos(qn(5)) 0 cos(qn(8)) 0 cos(qn(11)) 0 0 cos(qn(15)) 0 0];
    Sty = [0 0 sin(qn(5)) 0 sin(qn(8)) 0 sin(qn(11)) 0 0 cos(qn(15)) 0 0];
    Stz = [0 1 1 1 1 1 1 1 1 1 1 1];
    % Sy = [0 1 -1 -1 1 1 -1 -1 -1 1 1 1];
    % Sz = [0 0 0 0 0 0 0 0 0 0 0 0];
    H_joint = fk(qn,Ln);
    H_CoM = fk_CoM(fk(qn,Ln),Ln,bsip,80);
    v(:,1,i) = v_base(i,:);
    a(:,1,i) = a_base(i,:);
    k = 1;
    a_c_torque = zeros(30,78);
    for j = 2:12    % Forward pass for segment speeds and velocities
        Ri = H_CoM(1:3,1:3,j);
        Rli = H_CoM(1:3,1:3,lambda(j));
        R = Ri/Rli;
        r = Rli*(H_CoM(1:3,4,j)-H_CoM(1:3,4,lambda(j)));
        rx = cross([r r r],eye(3));
        X = [R, -R*rx; zeros(3,3), R];
        S = [0 0; 0 0; 0 0; 0 Sx(j); Sy(j) 0; Sz(j) 0];
        
        a_c_torque(2*k-1,(j-2)*6+1:(j-1)*6) = [0 0 0 0 0 Stz(j)];
        a_c_torque(2*k,(j-2)*6+1:(j-1)*6) = -[0 0 0 0 0 Stz(j)];
        k = k+1;
        if Stx(j) || Sty(j)
            a_c_torque(2*k-1,(j-2)*6+1:(j-1)*6) = [0 0 0 Stx(j) Sty(j) 0];
            a_c_torque(2*k,(j-2)*6+1:(j-1)*6) = -[0 0 0 Stx(j) Sty(j) 0];
            k = k+1;
        end
        v(:,j,i) = X*v(:,lambda(j),i) + S*qdjoint(:,j);
        a(:,j,i) = X*a(:,lambda(j),i) + S*qddjoint(:,j) - spatial_cross(v(:,j,i))'*(S*qdjoint(:,j));
    end
    for j = 1:12
        f((j-1)*6+1:(j-1)*6+3) = bsip.(fnames{j}).m*a(1:3,j,i);
        f((j-1)*6+4:(j-1)*6+6) = bsip.(fnames{j}).I*a(4:6,j,i) + cross(v(4:6,j,i),bsip.(fnames{j}).I*v(4:6,j,i));
        R = inv(H_CoM(1:3,1:3,j)); r = H_CoM(1:3,4,j);
        rx = cross([r r r],eye(3)); X = [R, zeros(3,3); rx*R, R];
        f((j-1)*6+1:(j-1)*6+6) = X*f((j-1)*6+1:(j-1)*6+6);
        f((j-1)*6+3) = f((j-1)*6+3) - 9.81*bsip.(fnames{j}).m;
    end
    A = zeros(72,78); R_all = zeros(78,78);
    for s = 1:12
        for j = 1:13
            if A_graph(s,j)
                r = H_joint(1:3,4,joint(j))-H_CoM(1:3,4,s);
                rx = cross([r r r],eye(3));
                A((s-1)*6+1:s*6,(j-1)*6+1:j*6) = A_graph(s,j)*[eye(3), zeros(3,3); rx, eye(3)];
            end
        end
    end
    for j = 1:11
        R = inv(H_joint(1:3,1:3,joint(j))); r = H_joint(1:3,4,joint(j));
        rx = cross([r r r],eye(3)); R_all((j-1)*6+1:j*6,(j-1)*6+1:j*6) = [R, zeros(3,3); rx*R, R];
    end
    for j = 12:13
        r = H_joint(1:3,4,joint(j)); rx = cross([r r r],eye(3));
        R_all((j-1)*6+1:j*6,(j-1)*6+1:j*6) = [eye(3), zeros(3,3); rx, eye(3)];
    end
    a_eq = A*R_all; b_eq = f;
    s = zeros(1,78); 
    s(:,69) = 1; s(:,75) = -1;
    alpha = 500; beta = 10^-5; gamma = 10^5; mi = 0.5;
    % proveriti b_c_torque
    b_c_torque = [143; 234; 92; 67; 67; 71; 77; 64; 67; 92; 71; 67; 46; 77; 185; 190; 190; 190; 168; 100; 126; 126; 190; 185; 190; 190; 100; 168; 126; 126];
    % a_c_torque = zeros(30,78);
    if (i < 125 || i > 186)
        H = 2*alpha*(s).'*s + beta*eye(78) + gamma*(a_eq).'*a_eq;
        g = (-2*gamma*b_eq.'*a_eq).';
        a_c_uni = zeros(2,78); a_c_uni(1,69) = -1; a_c_uni(2,75) = -1; b_c_uni = zeros(2,1); 
        a_c_fric = zeros(8,78); a_c_fric(1:4,67:69) = [1 0 -mi; -1 0 -mi; 0 1 -mi; 0 -1 -mi]; a_c_fric(5:8,73:75) = [1 0 -mi; -1 0 -mi; 0 1 -mi; 0 -1 -mi];
        b_c_fric = zeros(8,1); 
        % a_c = cat(1,a_c_uni,a_c_fric,a_c_torque); 
        % b_c = cat(1,b_c_uni,b_c_fric,b_c_torque);
        a_c = cat(1,a_c_uni,a_c_fric); 
        b_c = cat(1,b_c_uni,b_c_fric);
        fint(:,i+1) = quadprog(H,g,a_c,b_c,[],[],[],[],fint(:,i));
    else
        a_eq = a_eq(:,1:66); s = s(:,1:66);
        H = 2*alpha*(s).'*s + beta*eye(66) + gamma*(a_eq).'*a_eq;
        g = (-2*gamma*b_eq.'*a_eq).';
        a_c = a_c_torque(:,1:66);
        fint(1:66,i+1) = quadprog(H,g,[],[],[],[],[],[],fint(1:66,i));
        % fint(1:66,i+1) = quadprog(H,g,a_c,b_c_torque,[],[],[],[],fint(1:66,i));
    end
end

%%
close all
ind = input("ind = ");
figure();
plot(reshape(v(:,ind,:),[6,340]).')
legend('vx','vy','vz','wx','wy','wz')
ax = gca;
ax.ColorOrder = [1 0 0; 0 1 0; 0 0 1];
ax.LineStyleOrder = ["-","--"];

%%
ind = input("ind = ");
figure();
plot(reshape(a(:,ind,:),[6,340]).')
legend('ax','ay','az','alfax','alfay','alfaz')
ax = gca;
ax.ColorOrder = [1 0 0; 0 1 0; 0 0 1];
ax.LineStyleOrder = ["-","--"];
%%
ind = input("ind = ");
figure()
plot(fint((ind-1)*6+1:ind*6,:).')
legend('fx','fy','fz','nx','ny','nz')
ax = gca;
ax.ColorOrder = [1 0 0; 0 1 0; 0 0 1];
ax.LineStyleOrder = ["-","--"];
clear
close all
clc
    
%%

vid = [];
tic
lm = readmatrix('landmarks.csv');  % čitanje podataka

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
% 
bsip = HumanBSIP(80,Ln);

for i = 0:lm(end,1)
    q(i+1,:) = ik(lm(i*33+1:(i+1)*33,3:5),qn,Ln);
    qn = q(i+1,:).';
    scatter3(-lm(i*33+2:(i+1)*33,5), lm(i*33+2:(i+1)*33,3), -lm(i*33+2:(i+1)*33,4));
    hold on
    H_CoM = fk_CoM(fk(qn,Ln),qn,Ln,bsip,80);
    n = length(H_CoM);
    scatter3(reshape(H_CoM(1,4,1:n),[n,1]), ...
    reshape(H_CoM(2,4,1:n),[n,1]), ...
    reshape(H_CoM(3,4,1:n),[n,1]),'r')
    plotJoints(cat(3,fk(q(i+1,:),Ln),H_CoM));
    r_CoM(i+1,:) = H_CoM(1:3,4,13)-[-lm(i*33+32,5)-lm(i*33+33,5); lm(i*33+32,3)+lm(i*33+33,3); -lm(i*33+32,4)-lm(i*33+33,4)]/2;
    plot3(r_CoM(:,1),r_CoM(:,2),r_CoM(:,3))
    hold off
    fr = getframe(gcf);
    vid = [vid, fr];
end


for i = 1:18
    q(:,i) = smooth(q(:,i));
end

plotAngles(q)
toc

figure()
plot(-q(:,13))
hold on
plot(q(:,17))

figure()
plot(q(:,11))
hold on
plot(-q(:,15))

figure()
plot(-q(:,12))
hold on
plot(q(:,16))

t = linspace(0,1/90*size(q,1),size(q,1)).';
for i=1:18
 qd(:,i) = diff(q(:,i)) ./diff(t);
 qdd(:,i) = diff(qd(:,i)) ./ diff(t(1:end-1));
end
for i=1:3
    % v_CoM(:,i) = diff(lowpass(r_CoM(:,i),0.01)) ./ diff(t);
    v_CoM(:,i) = diff(r_CoM(:,i)) ./ diff(t);
end
v_CoM(:,4:6) = [qd(:,2), qd(:,3), qd(:,1)];
figure()
plot(t(1:end-1),v_CoM); legend(); hold on

%%
vz = [0; v_CoM(:,3)];
z_CoM = smooth(r_CoM(:,3));
[vz_max,i_max] = max(vz);
for i = i_max+1:181
    vz(i) = vz(i-1)-9.81/90;
    z_CoM(i) = z_CoM(i-1)+vz(i-1)/90;
end

figure()
plot(smooth(r_CoM(:,3))); hold on
plot(z_CoM)
plot(vz)

%%
close all
figure()
displacement = z_CoM - r_CoM(:,3);
for i = 0:lm(end,1)
    qn = q(i+1,:);
    disp = [(lm(i*33+32,5)+lm(i*33+33,5))/2;-(lm(i*33+32,3)+lm(i*33+33,3))/2;(lm(i*33+32,4)+lm(i*33+33,4))/2+displacement(i+1)];
    scatter3(-lm(i*33+2:(i+1)*33,5)+disp(1)*ones(32,1), ...
        lm(i*33+2:(i+1)*33,3)+disp(2)*ones(32,1), ...
        -lm(i*33+2:(i+1)*33,4)+disp(3)*ones(32,1));
    hold on
    H_CoM = fk_CoM(fk(qn,Ln),qn,Ln,bsip,80);
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
end

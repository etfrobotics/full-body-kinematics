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
Lh = (norm(lmT(16,:)-(lmT(18,:)+lmT(20,:))/2) + norm(lmT(15,:)-(lmT(17,:)+lmT(19,:))/2))/2;
Lbl = norm(lmT(23,:)-lmT(24,:))/2;
Ll1 = (norm(lmT(23,:)-lmT(25,:)) + norm(lmT(24,:)-lmT(26,:)))/2;
Ll2 = (norm(lmT(25,:)-lmT(27,:)) + norm(lmT(26,:)-lmT(28,:)))/2;
Lf1 = 0.8*(norm(lmT(27,:)-lmT(29,:)) + norm(lmT(28,:)-lmT(30,:)))/2;
Lf2 = 0.6*(norm(lmT(27,:)-lmT(29,:)) + norm(lmT(28,:)-lmT(30,:)))/2;
Lf3 = 0.6*(norm(lmT(29,:)-lmT(31,:)) + norm(lmT(30,:)-lmT(32,:)))/2;
Lf4 = 0.2*(norm(lmT(29,:)-lmT(31,:)) + norm(lmT(30,:)-lmT(32,:)))/2;

Ln = [Lbw;Lwt;Lta;La1;La2;Lh;Lbl;Ll1;Ll2;Lf1;Lf2;Lf3;Lf4];

for i = 0:lm(end,1)
    q(i+1,:) = ik(lm(i*33+1:(i+1)*33,3:5),Ln);
    scatter3(-lm(i*33+2:(i+1)*33,5), lm(i*33+2:(i+1)*33,3), -lm(i*33+2:(i+1)*33,4));
    hold on
    plotJoints(fk(q(i+1,:),Ln));
    hold off
    fr = getframe(gcf);
    vid = [vid, fr];
end

plotAngles(q)
toc
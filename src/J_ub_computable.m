function J_ub = J_ub_computable(in1,in2,in3)
%J_ub_computable
%    J_ub = J_ub_computable(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    06-Jul-2024 22:48:34

L1 = in2(1,:);
L2 = in2(2,:);
L3 = in2(3,:);
L4 = in2(4,:);
L5 = in2(5,:);
q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
q4 = in1(4,:);
q5 = in1(5,:);
q6 = in1(6,:);
q7 = in1(7,:);
q8 = in1(8,:);
xg1 = in3(1,:);
xg2 = in3(2,:);
xg3 = in3(3,:);
xg4 = in3(4,:);
xg5 = in3(5,:);
xg6 = in3(6,:);
xg7 = in3(7,:);
xg8 = in3(8,:);
xg9 = in3(9,:);
xg10 = in3(10,:);
xg11 = in3(11,:);
xg12 = in3(12,:);
xg13 = in3(13,:);
xg14 = in3(14,:);
xg15 = in3(15,:);
xg16 = in3(16,:);
xg17 = in3(17,:);
xg18 = in3(18,:);
t2 = cos(q1);
t3 = cos(q4);
t4 = cos(q5);
t5 = cos(q7);
t6 = cos(q8);
t7 = sin(q1);
t8 = sin(q4);
t9 = sin(q5);
t10 = sin(q7);
t11 = sin(q8);
t14 = pi./2.0;
t12 = L1.*t2;
t13 = L1.*t7;
t17 = -t14;
t18 = q6+t14;
t46 = t2.*6.123233995736766e-17;
t47 = t7.*6.123233995736766e-17;
t65 = t2.*3.749399456654644e-33;
t66 = t7.*3.749399456654644e-33;
t15 = -t12;
t16 = -t13;
t19 = q2+t17;
t20 = q3+t17;
t21 = cos(t18);
t22 = sin(t18);
t48 = -t47;
t67 = -t65;
t68 = -t66;
t23 = cos(t19);
t24 = cos(t20);
t25 = sin(t19);
t26 = sin(t20);
t27 = t2.*t23;
t28 = t2.*t25;
t29 = t7.*t23;
t30 = t7.*t25;
t49 = t23.*6.123233995736766e-17;
t50 = L2.*t25.*6.123233995736766e-17;
t58 = t22.*t25.*6.123233995736766e-17;
t59 = t25.*t26.*6.123233995736766e-17;
t69 = t23.*3.749399456654644e-33;
t79 = t21.*t25.*3.749399456654644e-33;
t80 = t24.*t25.*3.749399456654644e-33;
t93 = t23.*2.295845021658468e-49;
t31 = L2.*t27;
t32 = L2.*t28;
t33 = L2.*t29;
t34 = L2.*t30;
t35 = -t30;
t40 = t28+t29;
t51 = t27.*6.123233995736766e-17;
t52 = t28.*6.123233995736766e-17;
t53 = t29.*6.123233995736766e-17;
t54 = t30.*6.123233995736766e-17;
t70 = t27.*3.749399456654644e-33;
t71 = t28.*3.749399456654644e-33;
t72 = t29.*3.749399456654644e-33;
t73 = t30.*3.749399456654644e-33;
t74 = t69+1.0;
t76 = t49-6.123233995736766e-17;
t36 = -t31;
t37 = -t32;
t38 = -t33;
t39 = -t34;
t41 = t27+t35;
t42 = t22.*t40;
t43 = t26.*t40;
t55 = -t51;
t56 = -t52;
t57 = -t53;
t60 = t21.*t40.*6.123233995736766e-17;
t61 = t24.*t40.*6.123233995736766e-17;
t75 = -t73;
t77 = L3.*t74;
t82 = t21.*t76;
t83 = t24.*t76;
t85 = t40+t66;
t94 = t22.*t76.*6.123233995736766e-17;
t95 = t26.*t76.*6.123233995736766e-17;
t106 = t48+t52+t53;
t44 = t22.*t41;
t45 = t26.*t41;
t62 = -t61;
t63 = t21.*t41.*6.123233995736766e-17;
t64 = t24.*t41.*6.123233995736766e-17;
t78 = -t77;
t84 = -t82;
t86 = t41+t65;
t87 = t21.*t85;
t88 = t24.*t85;
t98 = -t95;
t105 = t46+t54+t55;
t109 = L3.*t106;
t110 = t22.*t85.*6.123233995736766e-17;
t111 = t26.*t85.*6.123233995736766e-17;
t115 = t59+t83;
t89 = t21.*t86;
t90 = t24.*t86;
t91 = -t88;
t97 = t44+t87;
t108 = L3.*t105;
t112 = -t110;
t113 = t22.*t86.*6.123233995736766e-17;
t114 = t26.*t86.*6.123233995736766e-17;
t116 = t58+t84;
t117 = L4.*t3.*t115;
t121 = t74+t80+t98;
t127 = t47+t56+t57+t64+t111;
t92 = -t89;
t99 = t43+t90;
t101 = L4.*t5.*t97;
t102 = t45+t91;
t118 = L4.*t5.*t116;
t123 = L4.*t8.*t121;
t124 = t60+t105+t113;
t125 = t63+t106+t112;
t126 = t62+t105+t114;
t130 = L4.*t8.*t127;
t100 = t42+t92;
t103 = L4.*t3.*t99;
t107 = L4.*t3.*t102;
t128 = L4.*t10.*t124;
t129 = L4.*t10.*t125;
t131 = L4.*t8.*t126;
t104 = L4.*t5.*t100;
t132 = -t131;
et1 = t93-t24.*t25.*6.123233995736766e-17+t26.*t76+t8.*t115.*6.123233995736766e-17-t3.*t121.*6.123233995736766e-17;
et2 = 6.123233995736766e-17;
et3 = t93+t21.*t25.*6.123233995736766e-17+t22.*t76-t10.*t116.*6.123233995736766e-17+t5.*(-t74+t79+t94).*6.123233995736766e-17;
et4 = 6.123233995736766e-17;
et5 = (t15+t32+t33+t108+xg3).^2+(t16+t34+t36+t109+xg10).^2+(-t50+t78+t118+xg14+L4.*t10.*(-t74+t79+t94)).^2+(t13+t31+t39+t107+t109+t130-xg4).^2+(t12+t37+t38+t104+t108+t128-xg15).^2+(t50+t78+t117+t123-xg8+L5.*t4.*(t3.*t115+t8.*t121)-L5.*t9.*(et1+et2)).^2;
et6 = (t15+t32+t33+t103+t108+t132+xg9+L5.*t4.*(t3.*t99-t8.*t126)-L5.*t9.*(t67+t70+t75-t24.*t40+t8.*t99.*6.123233995736766e-17+t26.*t86+t3.*t126.*6.123233995736766e-17)).^2;
et7 = (t16+t34+t36+t101+t109+t129+xg16+L5.*t6.*(t5.*t97+t10.*t125)-L5.*t11.*(t68+t71+t72-t21.*t41+t10.*t97.*6.123233995736766e-17+t22.*t85-t5.*t125.*6.123233995736766e-17)).^2+(-t50+t78+t118+xg17+L4.*t10.*(-t74+t79+t94)+L5.*t11.*(et3+et4)+L5.*t6.*(t5.*t116+t10.*(-t74+t79+t94))).^2+(t13+t31+t39+t109-xg1).^2;
et8 = (t15+t32+t33-t108+xg12).^2+(t50+t78+t117+t123-xg5).^2+(-t50+t77+xg2).^2+(t50+t77-xg11).^2+(t13+t31+t39+t107+t109+t130-xg7+L5.*t4.*(t3.*t102+t8.*t127)+L5.*t9.*(t68+t71+t72+t24.*t41-t8.*t102.*6.123233995736766e-17+t26.*t85+t3.*t127.*6.123233995736766e-17)).^2;
et9 = (t12+t37+t38+t104+t108+t128-xg18+L5.*t6.*(t5.*t100+t10.*t124)+L5.*t11.*(t67+t70+t75+t21.*t40+t22.*t86-t10.*t100.*6.123233995736766e-17+t5.*t124.*6.123233995736766e-17)).^2+(t15+t32+t33+t103+t108+t132+xg6).^2+(t16+t34+t36+t101+t109+t129+xg13).^2;
J_ub = et5+et6+et7+et8+et9;
end

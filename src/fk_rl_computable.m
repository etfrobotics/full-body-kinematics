function fk_rl = fk_rl_computable(in1,in2)
%FK_RL_COMPUTABLE
%    FK_RL = FK_RL_COMPUTABLE(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    22-Mar-2024 14:09:48

Ll1 = in2(1,:);
Ll2 = in2(2,:);
Ll3 = in2(3,:);
Ll4 = in2(4,:);
Ll5 = in2(5,:);
Ll6 = in2(6,:);
Ll7 = in2(7,:);
ql1 = in1(1,:);
ql2 = in1(2,:);
ql3 = in1(3,:);
ql4 = in1(4,:);
ql5 = in1(5,:);
ql6 = in1(6,:);
ql7 = in1(7,:);
t2 = cos(ql5);
t3 = cos(ql6);
t4 = sin(ql5);
t5 = sin(ql6);
t6 = Ll5+Ll6;
t7 = -Ll1;
t8 = pi./2.0;
t9 = -t8;
t10 = ql1+t8;
t11 = ql2+t8;
t12 = ql7+t8;
t13 = ql3+t9;
t14 = ql4+t9;
t15 = cos(t10);
t16 = cos(t11);
t17 = cos(t12);
t18 = sin(t10);
t19 = sin(t11);
t20 = sin(t12);
t21 = cos(t13);
t22 = cos(t14);
t23 = sin(t13);
t24 = sin(t14);
t25 = t15.*t16;
t26 = t15.*t19;
t27 = t16.*t18;
t28 = t18.*t19;
t33 = t15.*6.123233995736766e-17;
t34 = t16.*6.123233995736766e-17;
t35 = t18.*6.123233995736766e-17;
t56 = t15.*3.749399456654644e-33;
t57 = t16.*3.749399456654644e-33;
t58 = t16+3.749399456654644e-33;
t59 = t18.*3.749399456654644e-33;
t92 = t15.*2.295845021658468e-49;
t93 = t16.*2.295845021658468e-49;
t94 = t18.*2.295845021658468e-49;
t109 = t15.*1.405799628556214e-65;
t110 = t18.*1.405799628556214e-65;
t29 = t19.*t21;
t30 = t19.*t23;
t31 = -t26;
t36 = -t34;
t37 = -t35;
t38 = t25.*6.123233995736766e-17;
t39 = t26.*6.123233995736766e-17;
t40 = t27.*6.123233995736766e-17;
t41 = t28.*6.123233995736766e-17;
t60 = -t56;
t61 = -t57;
t62 = -t59;
t63 = Ll2.*t58;
t64 = t34-6.123233995736766e-17;
t66 = t25.*3.749399456654644e-33;
t67 = t26.*3.749399456654644e-33;
t68 = t27.*3.749399456654644e-33;
t69 = t28.*3.749399456654644e-33;
t95 = -t92;
t96 = -t94;
t97 = t25.*2.295845021658468e-49;
t98 = t26.*2.295845021658468e-49;
t99 = t27.*2.295845021658468e-49;
t100 = t28.*2.295845021658468e-49;
t117 = t25.*1.405799628556214e-65;
t118 = t27.*1.405799628556214e-65;
t32 = -t30;
t42 = -t38;
t43 = -t39;
t44 = -t40;
t45 = -t41;
t46 = t30.*6.123233995736766e-17;
t48 = t27+t39;
t65 = -t63;
t70 = -t67;
t71 = -t68;
t72 = -t69;
t73 = t30.*3.749399456654644e-33;
t74 = t21.*t64;
t75 = t23.*t64;
t79 = t26+t37+t40;
t101 = -t97;
t102 = -t99;
t105 = t18+t39+t68;
t47 = -t46;
t49 = t25+t45;
t50 = t21.*t48;
t51 = t23.*t48;
t76 = -t74;
t77 = t29+t75;
t78 = t28+t33+t42;
t84 = Ll2.*t79;
t103 = t74.*6.123233995736766e-17;
t106 = t15+t45+t66;
t111 = t21.*t105;
t112 = t23.*t105;
t119 = t74.*3.749399456654644e-33;
t52 = t21.*t49;
t53 = t23.*t49;
t54 = -t51;
t80 = t51.*6.123233995736766e-17;
t83 = Ll2.*t78;
t87 = -t84;
t88 = t22.*t77;
t89 = t24.*t77;
t104 = -t103;
t107 = t51.*3.749399456654644e-33;
t113 = t21.*t106;
t114 = t23.*t106;
t115 = -t111;
t116 = -t112;
t120 = t111.*6.123233995736766e-17;
t125 = t111.*3.749399456654644e-33;
t136 = t30+t36+t76-2.295845021658468e-49;
t55 = -t53;
t81 = t53.*6.123233995736766e-17;
t82 = -t80;
t86 = -t83;
t90 = Ll3.*t88;
t108 = t53.*3.749399456654644e-33;
t121 = t113.*6.123233995736766e-17;
t122 = -t120;
t124 = t50+t114;
t126 = t52+t116;
t127 = t113.*3.749399456654644e-33;
t135 = t46+t58+t104;
t162 = t45+t54+t60+t66+t113;
t85 = -t81;
t91 = -t90;
t123 = -t121;
t128 = t22.*t124;
t129 = t24.*t124;
t131 = t22.*t126;
t132 = t24.*t126;
t137 = t22.*t135;
t138 = t24.*t135;
t150 = t78+t82+t121;
t151 = t31+t35+t44+t81+t120;
t163 = t43+t55+t59+t71+t115;
t130 = Ll3.*t128;
t133 = Ll3.*t131;
t139 = Ll3.*t138;
t140 = -t138;
t141 = t89+t137;
t152 = t22.*t150;
t153 = t24.*t150;
t154 = t22.*t151;
t155 = t24.*t151;
t142 = t88+t140;
t143 = t2.*t141;
t144 = t4.*t141;
t156 = Ll3.*t153;
t157 = Ll3.*t155;
t158 = -t152;
t164 = t7+t65+t91+t139;
t165 = t128+t153;
t168 = t132+t154;
t178 = -t2.*(t131-t155);
t180 = -t4.*(t131-t155);
t181 = t4.*(t131-t155);
t145 = t2.*t142;
t146 = -t143;
t147 = t4.*t142;
t159 = -t157;
t160 = t143.*6.123233995736766e-17;
t166 = t129+t158;
t167 = t4.*t165;
t169 = t2.*t165;
t175 = t2.*t168;
t176 = t4.*t168;
t188 = t181.*(-6.123233995736766e-17);
t189 = t86+t130+t156;
t149 = -t147;
t161 = t147.*6.123233995736766e-17;
t170 = t4.*t166;
t171 = -t167;
t173 = t2.*t166;
t179 = -t175;
t182 = t167.*6.123233995736766e-17;
t186 = t175.*6.123233995736766e-17;
t190 = t87+t133+t159;
t192 = -t5.*(t144-t145);
t195 = t3.*(t144-t145);
t197 = t5.*(t144-t145).*(-6.123233995736766e-17);
t203 = t46+t61+t104+t143+t147-1.405799628556214e-65;
t211 = t176+t178;
t244 = t67+t81+t96+t99+t120+t175+t181;
t177 = -t173;
t183 = -t182;
t184 = t173.*6.123233995736766e-17;
t187 = -t186;
t196 = Ll4.*t195;
t198 = t32+t34+t74+t160+t161+2.295845021658468e-49;
t204 = t47+t57+t103+t146+t149+1.405799628556214e-65;
t210 = -t5.*(t169-t170);
t213 = t3.*(t169-t170);
t215 = t5.*t211;
t216 = t3.*t211;
t220 = t5.*(t169-t170).*(-6.123233995736766e-17);
t221 = t5.*(t169-t170).*6.123233995736766e-17;
t239 = t69+t80+t92+t101+t123+t167+t173;
t245 = Ll5.*t244;
t246 = t70+t85+t94+t102+t122+t179+t180;
t247 = t6.*t244;
t185 = -t184;
t199 = t3.*t198;
t200 = t5.*t198;
t205 = Ll5.*t204;
t206 = t6.*t204;
t214 = Ll4.*t213;
t217 = Ll4.*t216;
t218 = -t216;
t222 = t215.*6.123233995736766e-17;
t225 = t3.*(-t41-t51-t56+t66+t113+t182+t184);
t226 = t5.*(-t41-t51-t56+t66+t113+t182+t184);
t227 = t39+t53+t62+t68+t111+t187+t188;
t240 = Ll5.*t239;
t241 = t6.*t239;
t242 = t72+t82+t95+t97+t121+t171+t177;
t248 = -t247;
t201 = Ll4.*t200;
t202 = t199.*6.123233995736766e-17;
t207 = -t206;
t219 = -t217;
t223 = -t222;
t228 = Ll4.*t226;
t229 = -t225;
t230 = -t226;
t232 = t3.*t227;
t233 = t5.*t227;
t235 = t225.*6.123233995736766e-17;
t238 = t195+t200;
t243 = -t241;
t231 = -t228;
t234 = Ll4.*t233;
t236 = t232.*6.123233995736766e-17;
t250 = t213+t230;
t251 = t197+t202+t203;
t252 = t218+t233;
t253 = t221+t235+t239;
t237 = -t236;
t254 = t223+t237+t244;
et1 = t18;
et2 = 8.608040076769528e-82;
et3 = t27;
et4 = -8.608040076769528e-82;
et5 = et1.*et2-t26.*1.405799628556214e-65;
et6 = et3.*et4-t53.*2.295845021658468e-49;
et7 = t111.*(-2.295845021658468e-49);
et8 = t175.*(-3.749399456654644e-33);
et9 = t181.*(-3.749399456654644e-33)+t223+t237-t17.*t254+t20.*(t216-t233);
et10 = -t73+t93+t119-t160-t161+t192+t199-t17.*(-t202+t204+t5.*(t144-t145).*6.123233995736766e-17).*6.123233995736766e-17-t20.*t238.*6.123233995736766e-17;
et11 = 8.608040076769528e-82;
et12 = t73-t93-t119+t160+t161-t199+t5.*(t144-t145);
et13 = -8.608040076769528e-82;
et14 = t16.*(-1.405799628556214e-65);
et15 = t30.*2.295845021658468e-49;
et16 = t74.*(-2.295845021658468e-49);
et17 = t143.*3.749399456654644e-33+t147.*3.749399456654644e-33;
et18 = -t202-t17.*(-t202+t204+t5.*(t144-t145).*6.123233995736766e-17)-t20.*t238+t5.*(t144-t145).*6.123233995736766e-17;
et19 = -5.27090436347397e-98;
et20 = t15;
et21 = -8.608040076769528e-82;
et22 = t25;
et23 = 8.608040076769528e-82;
et24 = et20.*et21+et22.*et23;
et25 = t28.*(-1.405799628556214e-65);
et26 = t51.*(-2.295845021658468e-49);
et27 = t113.*2.295845021658468e-49;
et28 = t167.*(-3.749399456654644e-33);
et29 = t173.*(-3.749399456654644e-33)+t221+t235-t17.*t253-t20.*t250;
fk_rl = reshape([1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,-1.0,0.0,0.0,0.0,t7,0.0,1.0,t15,0.0,t18,0.0,t37,-1.0,t33,0.0,t18,-6.123233995736766e-17,-t15,0.0,0.0,t7,0.0,1.0,t49,-t19,t48,0.0,-t18+t43+t71,t36+6.123233995736766e-17,t106,0.0,t31+t35+t44,-t58,-t28-t33+t38,0.0,0.0,t7,0.0,1.0,t126,-t77,t124,0.0,t79+t85+t122,t135,t150,0.0,t163,t136,t162,0.0,t87,t7+t65,t86,1.0,t131-t155,-t88+t138,t165,0.0,-t168,t141,-t129+t152,0.0,t163,t136,t162,0.0,t190,t164,t189,1.0,-t176+t2.*(t131-t155),t144-t145,t169-t170,0.0,t227,t198,t41+t51+t56-t66-t113+t183+t185,0.0,t246,t203,t242,0.0,t190,t164,t189,1.0,t252,t238,t250,0.0,t215+t232,t192+t199,t210+t229,0.0,t246,t203,t242,0.0,t190+t219+t234+t245,t164+t196+t201+t205,t189+t214+t231+t240,1.0,t252,t238,t250,0.0,t222+t236+t246,t251,t220-t235+t242,0.0,-t98-t108+t110-t118-t125+t187+t188-t215-t232,et12+et13,-t100-t107-t109+t117+t127+t183+t185+t225+t5.*(t169-t170),0.0,t190+t219+t234+t245+t248,t164+t196+t201+t205+t207,t189+t214+t231+t240+t243,1.0,-t20.*t254-t17.*(t216-t233),-t20.*(-t202+t204+t5.*(t144-t145).*6.123233995736766e-17)+t17.*t238,t17.*t250-t20.*t253,0.0,t98+t108-t110+t118+t125+t181.*6.123233995736766e-17+t186+t215+t232-t17.*t254.*6.123233995736766e-17+t20.*(t216-t233).*6.123233995736766e-17,et10+et11,t100+t107+t109-t117-t127+t182+t184+t210+t229-t17.*t253.*6.123233995736766e-17-t20.*t250.*6.123233995736766e-17,0.0,et5+et6+et7+et8+et9,et14+et15+et16+et17+et18+et19,et24+et25+et26+et27+et28+et29,0.0,t190+t219+t234+t245+t248-Ll7.*t20.*t254-Ll7.*t17.*(t216-t233),t164+t196+t201+t205+t207-Ll7.*t20.*(-t202+t204+t5.*(t144-t145).*6.123233995736766e-17)+Ll7.*t17.*t238,t189+t214+t231+t240+t243+Ll7.*t17.*t250-Ll7.*t20.*t253,1.0],[4,4,9]);
end

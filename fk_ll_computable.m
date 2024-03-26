function fk_ll = fk_ll_computable(in1,in2)
%FK_LL_COMPUTABLE
%    FK_LL = FK_LL_COMPUTABLE(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    22-Mar-2024 14:09:53

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
t7 = pi./2.0;
t8 = -t7;
t9 = ql3+t7;
t10 = ql4+t7;
t11 = ql1+t8;
t12 = ql2+t8;
t13 = ql7+t8;
t14 = cos(t9);
t15 = cos(t10);
t16 = sin(t9);
t17 = sin(t10);
t18 = cos(t11);
t19 = cos(t12);
t20 = cos(t13);
t21 = sin(t11);
t22 = sin(t12);
t23 = sin(t13);
t24 = -t18;
t25 = -t21;
t26 = t14.*t22;
t27 = t16.*t22;
t28 = t18.*t19;
t29 = t18.*t22;
t30 = t19.*t21;
t31 = t21.*t22;
t35 = t18.*6.123233995736766e-17;
t36 = t19.*6.123233995736766e-17;
t37 = t21.*6.123233995736766e-17;
t57 = t18.*3.749399456654644e-33;
t58 = t19.*3.749399456654644e-33;
t59 = t19+3.749399456654644e-33;
t60 = t21.*3.749399456654644e-33;
t89 = t18.*2.295845021658468e-49;
t90 = t19.*2.295845021658468e-49;
t91 = t21.*2.295845021658468e-49;
t107 = t18.*1.405799628556214e-65;
t108 = t21.*1.405799628556214e-65;
t32 = -t27;
t33 = t22.*t24;
t34 = t22.*t25;
t38 = -t35;
t39 = -t37;
t40 = t27.*6.123233995736766e-17;
t41 = t28.*6.123233995736766e-17;
t42 = t29.*6.123233995736766e-17;
t43 = t30.*6.123233995736766e-17;
t44 = t31.*6.123233995736766e-17;
t61 = -t57;
t62 = -t58;
t63 = -t60;
t64 = Ll2.*t59;
t65 = t36-6.123233995736766e-17;
t66 = t27.*3.749399456654644e-33;
t67 = t28.*3.749399456654644e-33;
t68 = t29.*3.749399456654644e-33;
t69 = t30.*3.749399456654644e-33;
t70 = t31.*3.749399456654644e-33;
t92 = -t89;
t93 = -t91;
t94 = t28.*2.295845021658468e-49;
t95 = t29.*2.295845021658468e-49;
t96 = t30.*2.295845021658468e-49;
t97 = t31.*2.295845021658468e-49;
t109 = -t107;
t110 = -t108;
t116 = t28.*1.405799628556214e-65;
t117 = t30.*1.405799628556214e-65;
t45 = -t40;
t46 = -t41;
t47 = -t42;
t48 = -t43;
t49 = -t44;
t50 = t30+t42;
t71 = -t66;
t72 = -t70;
t73 = t14.*t65;
t74 = t16.*t65;
t83 = t29+t39+t43;
t98 = -t94;
t99 = -t97;
t102 = t21+t42+t69;
t51 = t28+t49;
t52 = t14.*t50;
t53 = t16.*t50;
t75 = t26+t74;
t81 = t31+t35+t46;
t86 = Ll2.*t83;
t100 = t73.*6.123233995736766e-17;
t104 = t18+t49+t67;
t111 = t14.*t102;
t112 = t16.*t102;
t118 = t73.*3.749399456654644e-33;
t135 = t32+t36+t73+2.295845021658468e-49;
t54 = t14.*t51;
t55 = t16.*t51;
t56 = -t53;
t76 = t53.*6.123233995736766e-17;
t77 = t15.*t75;
t78 = t17.*t75;
t85 = Ll2.*t81;
t101 = -t100;
t103 = t53.*3.749399456654644e-33;
t113 = t14.*t104;
t114 = t16.*t104;
t115 = -t112;
t119 = t111.*6.123233995736766e-17;
t123 = t111.*3.749399456654644e-33;
t79 = t55.*6.123233995736766e-17;
t80 = -t76;
t82 = Ll3.*t77;
t88 = -t85;
t105 = t55.*3.749399456654644e-33;
t106 = -t103;
t120 = t113.*6.123233995736766e-17;
t121 = -t119;
t124 = t113.*3.749399456654644e-33;
t125 = t52+t114;
t126 = t54+t115;
t136 = t40+t59+t101;
t165 = t42+t55+t63+t69+t111;
t166 = t49+t56+t61+t67+t113;
t84 = -t79;
t87 = -t82;
t122 = -t120;
t127 = t15.*t125;
t128 = t17.*t125;
t130 = t15.*t126;
t131 = t17.*t126;
t137 = t15.*t136;
t138 = t17.*t136;
t150 = t80+t81+t120;
t151 = t33+t37+t48+t79+t119;
t129 = Ll3.*t127;
t133 = Ll3.*t130;
t139 = Ll3.*t138;
t140 = -t138;
t141 = t78+t137;
t156 = t15.*t150;
t157 = t17.*t150;
t158 = t15.*t151;
t159 = t17.*t151;
t132 = -t129;
t142 = t77+t140;
t143 = t2.*t141;
t144 = t4.*t141;
t160 = Ll3.*t157;
t161 = Ll3.*t159;
t162 = -t156;
t167 = Ll1+t64+t87+t139;
t168 = t127+t157;
t171 = t131+t158;
t181 = -t2.*(t130-t159);
t182 = -t4.*(t130-t159);
t183 = t4.*(t130-t159);
t145 = t2.*t142;
t146 = -t143;
t147 = t4.*t142;
t152 = t143.*6.123233995736766e-17;
t163 = -t160;
t164 = -t161;
t169 = t128+t162;
t170 = t4.*t168;
t172 = t2.*t168;
t178 = t2.*t171;
t179 = t4.*t171;
t190 = t183.*(-6.123233995736766e-17);
t191 = t183.*6.123233995736766e-17;
t149 = -t147;
t153 = -t152;
t154 = t147.*6.123233995736766e-17;
t173 = t4.*t169;
t174 = -t170;
t176 = t2.*t169;
t177 = -t172;
t184 = t170.*6.123233995736766e-17;
t188 = t178.*6.123233995736766e-17;
t193 = -t5.*(t144-t145);
t196 = t3.*(t144-t145);
t198 = t86+t133+t164;
t199 = t88+t132+t163;
t200 = t5.*(t144-t145).*(-6.123233995736766e-17);
t213 = t179+t181;
t242 = t68+t79+t93+t96+t119+t178+t183;
t155 = -t154;
t180 = -t176;
t185 = -t184;
t186 = t176.*6.123233995736766e-17;
t189 = -t188;
t197 = Ll4.*t196;
t201 = t135+t152+t154;
t206 = t45+t58+t100+t146+t149+1.405799628556214e-65;
t210 = t173+t177;
t211 = -t3.*(t172-t173);
t215 = t5.*(t172-t173);
t216 = t5.*t213;
t217 = t3.*t213;
t230 = -t3.*(-t165+t188+t191);
t231 = -t5.*(-t165+t188+t191);
t235 = t3.*(-t165+t188+t191).*(-6.123233995736766e-17);
t236 = t3.*(-t165+t188+t191).*6.123233995736766e-17;
t237 = t70+t76+t89+t98+t122+t170+t176;
t243 = Ll5.*t242;
t245 = t6.*t242;
t187 = -t186;
t202 = t3.*t201;
t203 = t5.*t201;
t207 = Ll5.*t206;
t208 = t6.*t206;
t214 = Ll4.*t211;
t218 = Ll4.*t217;
t219 = -t217;
t222 = t215.*6.123233995736766e-17;
t223 = t216.*6.123233995736766e-17;
t225 = t166+t184+t186;
t226 = t3.*(-t44-t53-t57+t67+t113+t184+t186);
t227 = t5.*(-t44-t53-t57+t67+t113+t184+t186);
t228 = t165+t189+t190;
t232 = Ll4.*t231;
t238 = Ll5.*t237;
t239 = t6.*t237;
t240 = t72+t80+t92+t94+t120+t174+t180;
t244 = -t243;
t204 = Ll4.*t203;
t205 = t202.*6.123233995736766e-17;
t209 = -t207;
t220 = -t218;
t224 = -t223;
t229 = Ll4.*t227;
t233 = t196+t203;
t234 = t226.*6.123233995736766e-17;
t241 = -t239;
t247 = t211+t227;
t248 = t219+t231;
t246 = t40+t62+t101+t143+t147+t200+t205-1.405799628556214e-65;
t249 = t222+t234+t237;
t250 = t224+t236+t242;
et1 = t21;
et2 = -8.608040076769528e-82;
et3 = t30;
et4 = 8.608040076769528e-82;
et5 = et1.*et2+t29.*1.405799628556214e-65+et3.*et4;
et6 = t55.*2.295845021658468e-49;
et7 = t111.*2.295845021658468e-49;
et8 = t178.*3.749399456654644e-33+t183.*3.749399456654644e-33;
et9 = t223+t235+t20.*t250-t23.*(t217+t5.*(-t165+t188+t191));
et10 = t71+t90+t118+t153+t155+t193+t202-t20.*(-t205+t206+t5.*(t144-t145).*6.123233995736766e-17).*6.123233995736766e-17-t23.*t233.*6.123233995736766e-17;
et11 = 8.608040076769528e-82;
et12 = t71+t90+t118+t153+t155+t193+t202;
et13 = 8.608040076769528e-82;
et14 = t19.*1.405799628556214e-65;
et15 = t27.*(-2.295845021658468e-49);
et16 = t73.*2.295845021658468e-49;
et17 = t143.*(-3.749399456654644e-33);
et18 = t147.*(-3.749399456654644e-33)+t200+t205+t20.*(-t205+t206+t5.*(t144-t145).*6.123233995736766e-17)+t23.*t233;
et19 = 5.27090436347397e-98;
et20 = t18;
et21 = -8.608040076769528e-82;
et22 = t28;
et23 = 8.608040076769528e-82;
et24 = et20.*et21+et22.*et23;
et25 = t31.*(-1.405799628556214e-65);
et26 = t53.*(-2.295845021658468e-49);
et27 = t113.*2.295845021658468e-49;
et28 = t170.*(-3.749399456654644e-33);
et29 = t176.*(-3.749399456654644e-33)+t222+t234+t23.*(t227-t3.*(t172-t173))-t20.*t249;
fk_ll = reshape([1.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,Ll1,0.0,1.0,t18,0.0,t25,0.0,t39,-1.0,t38,0.0,t25,6.123233995736766e-17,t24,0.0,0.0,Ll1,0.0,1.0,t51,-t22,t47+t19.*t25,0.0,t25+t47-t69,-t65,t24+t44-t67,0.0,t83,t59,t34+t38+t41,0.0,0.0,Ll1,0.0,1.0,t126,-t75,-t125,0.0,t83+t84+t121,t136,t34+t38+t41+t76+t122,0.0,t165,t135,t166,0.0,t86,Ll1+t64,t88,1.0,t130-t159,-t77+t138,-t168,0.0,-t171,t141,t169,0.0,t165,t135,t166,0.0,t198,t167,t199,1.0,-t179+t2.*(t130-t159),t144-t145,t210,0.0,t228,t201,t225,0.0,t242,t206,t240,0.0,t198,t167,t199,1.0,t248,t233,t247,0.0,t216+t230,t193+t202,t215+t226,0.0,t242,t206,t240,0.0,t198+t220+t232+t244,t167+t197+t204+t209,t199+t214+t229+t238,1.0,t248,t233,t247,0.0,-t68+t84+t91-t96+t121-t178+t182+t223+t235,t246,t249,0.0,t95+t105+t110+t117+t123+t188+t191+t216+t230,et12+et13,t99+t106+t109+t116+t124+t185+t187+t215+t226,0.0,t198+t220+t232+t244+t245,t167+t197+t204+t208+t209,t199+t214+t229+t238+t241,1.0,-t23.*t250-t20.*(t217+t5.*(-t165+t188+t191)),-t23.*(-t205+t206+t5.*(t144-t145).*6.123233995736766e-17)+t20.*t233,t20.*(t227-t3.*(t172-t173))+t23.*t249,0.0,t95+t105+t110+t117+t123+t188+t191+t216+t230-t20.*t250.*6.123233995736766e-17+t23.*(t217+t5.*(-t165+t188+t191)).*6.123233995736766e-17,et10+et11,t99+t106+t109+t116+t124+t185+t187+t215+t226-t23.*(t227-t3.*(t172-t173)).*6.123233995736766e-17+t20.*t249.*6.123233995736766e-17,0.0,et5+et6+et7+et8+et9,et14+et15+et16+et17+et18+et19,et24+et25+et26+et27+et28+et29,0.0,t198+t220+t232+t244+t245-Ll7.*t23.*t250-Ll7.*t20.*(t217+t5.*(-t165+t188+t191)),t167+t197+t204+t208+t209-Ll7.*t23.*(-t205+t206+t5.*(t144-t145).*6.123233995736766e-17)+Ll7.*t20.*t233,t199+t214+t229+t238+t241+Ll7.*t20.*(t227-t3.*(t172-t173))+Ll7.*t23.*t249,1.0],[4,4,9]);
end
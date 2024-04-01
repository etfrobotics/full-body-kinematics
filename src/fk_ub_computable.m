function fk_ub = fk_ub_computable(in1,in2)
%FK_UB_COMPUTABLE
%    FK_UB = FK_UB_COMPUTABLE(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    20-Mar-2024 09:16:51

L1 = in2(1,:);
L2 = in2(2,:);
L3 = in2(3,:);
L4 = in2(4,:);
L5 = in2(5,:);
L6 = in2(6,:);
q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
q4 = in1(4,:);
q5 = in1(5,:);
q6 = in1(6,:);
q7 = in1(7,:);
q8 = in1(8,:);
q9 = in1(9,:);
q10 = in1(10,:);
q11 = in1(11,:);
q12 = in1(12,:);
q13 = in1(13,:);
q14 = in1(14,:);
q15 = in1(15,:);
q16 = in1(16,:);
q17 = in1(17,:);
q18 = in1(18,:);
t2 = cos(q1);
t3 = cos(q4);
t4 = cos(q8);
t5 = cos(q11);
t6 = cos(q15);
t7 = cos(q18);
t8 = sin(q1);
t9 = sin(q4);
t10 = sin(q8);
t11 = sin(q11);
t12 = sin(q15);
t13 = sin(q18);
t17 = pi./2.0;
t14 = L1.*t2;
t15 = L1.*t8;
t16 = -t8;
t18 = -t17;
t19 = q2+t17;
t20 = q3+t17;
t21 = q5+t17;
t22 = q6+t17;
t23 = q9+t17;
t24 = q10+t17;
t25 = q14+t17;
t56 = t2.*6.123233995736766e-17;
t57 = t8.*6.123233995736766e-17;
t81 = t2.*3.749399456654644e-33;
t82 = t8.*3.749399456654644e-33;
t103 = t2.*2.295845021658468e-49;
t104 = t8.*2.295845021658468e-49;
t26 = q7+t18;
t27 = q12+t18;
t28 = q13+t18;
t29 = q16+t18;
t30 = q17+t18;
t31 = cos(t19);
t32 = cos(t20);
t33 = cos(t21);
t34 = cos(t22);
t35 = cos(t23);
t36 = cos(t24);
t37 = cos(t25);
t38 = sin(t19);
t39 = sin(t20);
t40 = sin(t21);
t41 = sin(t22);
t42 = sin(t23);
t43 = sin(t24);
t44 = sin(t25);
t58 = -t57;
t59 = t2+t57;
t60 = t16+t56;
t83 = -t82;
t105 = -t103;
t45 = cos(t26);
t46 = cos(t27);
t47 = cos(t28);
t48 = cos(t29);
t49 = cos(t30);
t50 = sin(t26);
t51 = sin(t27);
t52 = sin(t28);
t53 = sin(t29);
t54 = sin(t30);
t55 = -t38;
t61 = t38.*6.123233995736766e-17;
t63 = t31.*t56;
t64 = t38.*t56;
t65 = t31.*t57;
t67 = t8.*t38.*(-6.123233995736766e-17);
t68 = t31.*t59;
t69 = t38.*t59;
t72 = t31.*t60;
t75 = t31.*3.749399456654644e-33;
t76 = t38.*3.749399456654644e-33;
t84 = t38.*t81;
t85 = t38.*t82;
t86 = t2.*t38.*(-3.749399456654644e-33);
t87 = t8.*t38.*(-3.749399456654644e-33);
t92 = t31.*2.295845021658468e-49;
t62 = -t61;
t70 = -t68;
t71 = t55.*t59;
t74 = t55.*t60;
t77 = -t75;
t78 = t31+t76;
t88 = t68.*6.123233995736766e-17;
t90 = t72.*6.123233995736766e-17;
t93 = -t92;
t79 = t32.*t78;
t80 = t39.*t78;
t89 = -t88;
t91 = -t90;
t94 = t65+t71;
t95 = t63+t74;
t109 = t61+t93+6.123233995736766e-17;
t113 = t60+t85+t88;
t118 = -t32.*(-t59+t84+t90);
t119 = -t39.*(-t59+t84+t90);
t121 = t32.*(-t59+t84+t90);
t96 = t39.*t94;
t97 = t32.*t94;
t100 = t32.*t95;
t101 = t39.*t95;
t102 = t80.*6.123233995736766e-17;
t110 = t32.*t109;
t111 = t39.*t109;
t114 = t59+t86+t91;
t116 = t39.*t113;
t117 = t32.*t113;
t128 = t121.*(-6.123233995736766e-17);
t98 = -t96;
t99 = -t97;
t106 = t96.*6.123233995736766e-17;
t107 = t101.*6.123233995736766e-17;
t112 = -t111;
t115 = t110.*6.123233995736766e-17;
t120 = -t117;
t127 = t117.*6.123233995736766e-17;
t132 = -t3.*(t97-t116);
t133 = t100+t119;
t134 = -t9.*(t97-t116);
t137 = t9.*(t97-t116);
t168 = t61+t80+t93+t110-2.295845021658468e-49;
t219 = t82+t85+t88+t96+t105+t117;
t231 = t81+t84+t90+t101+t104+t121;
t108 = -t107;
t122 = t79+t112;
t131 = t99+t116;
t138 = t3.*t133;
t139 = t9.*t133;
t143 = t137.*(-6.123233995736766e-17);
t144 = t137.*6.123233995736766e-17;
t149 = t137.*(-3.749399456654644e-33);
t150 = t137.*3.749399456654644e-33;
t155 = t137.*(-2.295845021658468e-49);
t157 = t55+t75+t102+t115+3.749399456654644e-33;
t163 = t137.*(-1.405799628556214e-65);
t165 = t137.*1.405799628556214e-65;
t169 = L2.*t168;
t172 = t33.*t168;
t173 = t40.*t168;
t174 = t46.*t168;
t175 = t51.*t168;
et1 = t137;
et2 = -8.608040076769528e-82;
t178 = et1.*et2;
et3 = t137;
et4 = 8.608040076769528e-82;
t179 = et3.*et4;
et5 = t137;
et6 = -5.27090436347397e-98;
t188 = et5.*et6;
et7 = t137;
et8 = 5.27090436347397e-98;
t189 = et7.*et8;
t199 = t58+t67+t70+t81+t106+t127;
t220 = L2.*t219;
t225 = t83+t87+t89+t98+t103+t120;
t229 = t33.*t219;
t230 = t40.*t219;
t233 = t46.*t219;
t234 = t51.*t219;
t235 = L2.*t231;
t240 = t33.*t231;
t241 = t40.*t231;
t242 = t46.*t231;
t243 = t51.*t231;
t123 = t3.*t122;
t124 = t9.*t122;
t140 = -t139;
t145 = t139.*6.123233995736766e-17;
t151 = t139.*3.749399456654644e-33;
t156 = t139.*2.295845021658468e-49;
t158 = t3.*t157;
t159 = t9.*t157;
t166 = t139.*1.405799628556214e-65;
t176 = -t172;
t177 = -t174;
et9 = t139;
et10 = 8.608040076769528e-82;
t180 = et9.*et10;
t184 = t172.*6.123233995736766e-17;
t186 = t174.*6.123233995736766e-17;
et11 = t139;
et12 = 5.27090436347397e-98;
t190 = et11.*et12;
t194 = t172.*3.749399456654644e-33;
t196 = t174.*3.749399456654644e-33;
t200 = t172.*2.295845021658468e-49;
t202 = t174.*2.295845021658468e-49;
t203 = t3.*t199;
t204 = t9.*t199;
t206 = t56+t64+t72+t82+t108+t128;
t212 = t172.*1.405799628556214e-65;
t214 = t174.*1.405799628556214e-65;
t223 = -t220;
et13 = t172;
et14 = 8.608040076769528e-82;
t224 = et13.*et14;
et15 = t174;
et16 = 8.608040076769528e-82;
t227 = et15.*et16;
t236 = -t230;
t244 = -t242;
t252 = t229.*6.123233995736766e-17;
t255 = t233.*6.123233995736766e-17;
t265 = t240.*6.123233995736766e-17;
t267 = t242.*6.123233995736766e-17;
t271 = t229.*3.749399456654644e-33;
t272 = t233.*3.749399456654644e-33;
t277 = t240.*3.749399456654644e-33;
t281 = t242.*3.749399456654644e-33;
t285 = t229.*2.295845021658468e-49;
t287 = t233.*2.295845021658468e-49;
t291 = t240.*2.295845021658468e-49;
t295 = t242.*2.295845021658468e-49;
t299 = t229.*1.405799628556214e-65;
t301 = t233.*1.405799628556214e-65;
t304 = t240.*1.405799628556214e-65;
t306 = t242.*1.405799628556214e-65;
et17 = t229;
et18 = 8.608040076769528e-82;
t311 = et17.*et18;
et19 = t233;
et20 = 8.608040076769528e-82;
t312 = et19.*et20;
et21 = t240;
et22 = 8.608040076769528e-82;
t314 = et21.*et22;
et23 = t242;
et24 = 8.608040076769528e-82;
t318 = et23.*et24;
t126 = -t124;
t129 = t124.*6.123233995736766e-17;
t135 = t124.*3.749399456654644e-33;
t141 = t124.*2.295845021658468e-49;
t146 = -t145;
t147 = t124.*1.405799628556214e-65;
t152 = -t151;
et25 = t124;
et26 = 8.608040076769528e-82;
t153 = et25.*et26;
t160 = -t158;
t161 = -t159;
et27 = t124;
et28 = 5.27090436347397e-98;
t162 = et27.*et28;
t167 = -t166;
t170 = t158.*6.123233995736766e-17;
t181 = -t180;
t182 = t158.*3.749399456654644e-33;
t185 = -t184;
t187 = -t186;
t191 = -t190;
t192 = t158.*2.295845021658468e-49;
t195 = -t194;
t197 = t158.*1.405799628556214e-65;
t201 = -t200;
t205 = -t203;
t207 = t3.*t206;
t208 = t9.*t206;
et29 = t158;
et30 = 8.608040076769528e-82;
t210 = et29.*et30;
t213 = -t212;
t215 = t203.*6.123233995736766e-17;
et31 = t158;
et32 = 5.27090436347397e-98;
t221 = et31.*et32;
t226 = -t224;
t228 = -t227;
t232 = t203.*3.749399456654644e-33;
t245 = t124+t158;
t251 = t203.*2.295845021658468e-49;
t254 = -t252;
t256 = -t255;
t268 = -t267;
t269 = t203.*1.405799628556214e-65;
t274 = -t272;
t280 = -t277;
et33 = t203;
et34 = 8.608040076769528e-82;
t283 = et33.*et34;
t286 = -t285;
t294 = -t291;
et35 = t203;
et36 = 5.27090436347397e-98;
t297 = et35.*et36;
t300 = -t299;
t305 = -t304;
t317 = -t314;
t320 = -t318;
t324 = t132+t204;
t325 = t137+t203;
t329 = t33.*(t204-t3.*(t97-t116));
t330 = t40.*(t204-t3.*(t97-t116));
t332 = t46.*(t204-t3.*(t97-t116));
t333 = t51.*(t204-t3.*(t97-t116));
t130 = -t129;
t136 = -t135;
t142 = -t141;
t148 = -t147;
t154 = -t153;
t164 = -t162;
t171 = -t170;
t183 = -t182;
t193 = -t192;
t198 = -t197;
t209 = -t207;
t211 = -t210;
t216 = -t215;
t217 = t207.*6.123233995736766e-17;
t222 = -t221;
t238 = t207.*3.749399456654644e-33;
t246 = L3.*t245;
t248 = t123+t161;
t250 = t126+t160;
t253 = -t251;
t257 = t207.*2.295845021658468e-49;
t270 = -t269;
t273 = t207.*1.405799628556214e-65;
t284 = -t283;
et37 = t207;
et38 = 8.608040076769528e-82;
t288 = et37.*et38;
t298 = -t297;
et39 = t207;
et40 = 5.27090436347397e-98;
t302 = et39.*et40;
t326 = L3.*t325;
t328 = t138+t208;
t331 = t140+t207;
t334 = -t333;
t344 = t330.*6.123233995736766e-17;
t346 = t333.*6.123233995736766e-17;
t352 = t330.*3.749399456654644e-33;
t353 = t333.*3.749399456654644e-33;
t356 = t330.*2.295845021658468e-49;
t358 = t333.*2.295845021658468e-49;
t363 = t330.*1.405799628556214e-65;
t365 = t333.*1.405799628556214e-65;
et41 = t330;
et42 = 8.608040076769528e-82;
t370 = et41.*et42;
et43 = t333;
et44 = 8.608040076769528e-82;
t371 = et43.*et44;
t403 = t236+t329;
t404 = t234+t332;
t407 = -t34.*(t230-t329);
t408 = -t41.*(t230-t329);
t411 = t34.*(t230-t329);
t420 = t41.*(t230-t329).*(-6.123233995736766e-17);
t421 = t41.*(t230-t329).*6.123233995736766e-17;
t428 = t41.*(t230-t329).*(-3.749399456654644e-33);
t429 = t41.*(t230-t329).*3.749399456654644e-33;
t434 = t41.*(t230-t329).*(-2.295845021658468e-49);
t435 = t41.*(t230-t329).*2.295845021658468e-49;
t441 = t41.*(t230-t329).*1.405799628556214e-65;
t218 = -t217;
t239 = -t238;
t249 = -t246;
t258 = -t257;
t260 = -t40.*t248;
t261 = -t46.*t248;
t262 = -t51.*t248;
t263 = t33.*t248;
t264 = t40.*t248;
t266 = t51.*t248;
t275 = -t273;
t289 = -t288;
t303 = -t302;
t327 = -t326;
t335 = -L3.*(t139+t209);
t336 = L3.*(t139+t209);
t337 = t33.*t328;
t338 = t40.*t328;
t339 = t46.*t328;
t340 = t51.*t328;
t345 = -t344;
t347 = -t346;
t357 = -t356;
t359 = -t358;
t364 = -t363;
t366 = -t365;
t372 = -t371;
t376 = t169+t246;
t401 = t14+t223+t326;
t409 = t47.*t404;
t410 = t52.*t404;
t467 = t252+t325+t344;
t468 = t256+t325+t346;
t276 = t264.*(-6.123233995736766e-17);
t278 = t264.*6.123233995736766e-17;
t279 = t266.*(-6.123233995736766e-17);
t282 = t266.*6.123233995736766e-17;
t290 = t264.*(-3.749399456654644e-33);
t292 = t264.*3.749399456654644e-33;
t293 = t266.*(-3.749399456654644e-33);
t296 = t266.*3.749399456654644e-33;
t308 = t264.*2.295845021658468e-49;
t309 = t266.*(-2.295845021658468e-49);
t310 = t266.*2.295845021658468e-49;
t315 = t264.*1.405799628556214e-65;
t316 = t266.*(-1.405799628556214e-65);
t319 = t266.*1.405799628556214e-65;
et45 = t264;
et46 = 8.608040076769528e-82;
t322 = et45.*et46;
et47 = t266;
et48 = -8.608040076769528e-82;
t323 = et47.*et48;
t341 = -t338;
t342 = -t339;
t343 = -t340;
t348 = t338.*6.123233995736766e-17;
t350 = t340.*6.123233995736766e-17;
t354 = t338.*3.749399456654644e-33;
t355 = t340.*3.749399456654644e-33;
t360 = t338.*2.295845021658468e-49;
t362 = t340.*2.295845021658468e-49;
t367 = t338.*1.405799628556214e-65;
t369 = t340.*1.405799628556214e-65;
et49 = t338;
et50 = 8.608040076769528e-82;
t373 = et49.*et50;
et51 = t340;
et52 = 8.608040076769528e-82;
t374 = et51.*et52;
t377 = t169+t249;
t378 = t173+t263;
t379 = t175+t261;
t402 = t14+t223+t327;
t405 = t15+t235+t335;
t406 = t15+t235+t336;
t413 = t241+t337;
t422 = t410.*6.123233995736766e-17;
t430 = t410.*3.749399456654644e-33;
t436 = t410.*2.295845021658468e-49;
t442 = t410.*1.405799628556214e-65;
t469 = t34.*t467;
t470 = t41.*t467;
t471 = t47.*t468;
t472 = t52.*t468;
t349 = -t348;
t351 = -t350;
t361 = -t360;
t375 = -t374;
t380 = t34.*t378;
t381 = t41.*t378;
t382 = t47.*t379;
t383 = t52.*t379;
t414 = t243+t342;
t415 = t34.*t413;
t416 = t41.*t413;
t423 = -t422;
t447 = t184+t245+t276;
t448 = t187+t245+t279;
t474 = t267+t331+t350;
t481 = t469.*6.123233995736766e-17;
t483 = t471.*6.123233995736766e-17;
t488 = t469.*3.749399456654644e-33;
t490 = t471.*3.749399456654644e-33;
t495 = t469.*2.295845021658468e-49;
t497 = t471.*2.295845021658468e-49;
t502 = t469.*1.405799628556214e-65;
t504 = t471.*1.405799628556214e-65;
t528 = t411+t470;
t533 = -t44.*(t409-t472);
t543 = t44.*(t409-t472).*(-6.123233995736766e-17);
t544 = t44.*(t409-t472).*6.123233995736766e-17;
t551 = t44.*(t409-t472).*(-3.749399456654644e-33);
t552 = t44.*(t409-t472).*3.749399456654644e-33;
t558 = t44.*(t409-t472).*2.295845021658468e-49;
t608 = t150+t232+t254+t345+t408+t469;
fk_ub = ft_1({L4,L5,L6,t10,t101,t102,t104,t11,t110,t113,t114,t115,t118,t12,t122,t129,t13,t130,t131,t133,t134,t135,t136,t137,t139,t14,t141,t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152,t153,t154,t155,t156,t16,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t195,t196,t197,t198,t199,t2,t201,t202,t205,t206,t209,t210,t211,t213,t214,t215,t216,t217,t218,t219,t221,t222,t223,t225,t226,t228,t229,t230,t231,t232,t233,t235,t238,t239,t240,t242,t243,t244,t245,t248,t250,t251,t252,t253,t254,t255,t257,t258,t260,t262,t264,t265,t266,t267,t268,t269,t270,t271,t272,t273,t274,t275,t276,t277,t278,t279,t280,t281,t282,t283,t284,t285,t286,t287,t288,t289,t290,t291,t292,t293,t294,t295,t296,t297,t298,t299,t300,t301,t302,t303,t304,t305,t306,t308,t309,t31,t310,t311,t312,t315,t316,t317,t319,t320,t322,t323,t324,t325,t328,t329,t330,t331,t334,t338,t339,t34,t340,t341,t343,t344,t345,t347,t348,t349,t35,t350,t351,t352,t353,t354,t355,t356,t357,t358,t359,t36,t360,t361,t362,t363,t364,t365,t366,t367,t369,t37,t370,t372,t373,t375,t376,t377,t378,t38,t380,t381,t382,t383,t4,t401,t402,t403,t404,t405,t406,t407,t409,t41,t410,t413,t414,t415,t416,t42,t420,t421,t422,t423,t428,t429,t43,t430,t434,t435,t436,t44,t441,t442,t447,t448,t45,t46,t469,t47,t470,t471,t472,t474,t48,t481,t483,t488,t49,t490,t495,t497,t5,t50,t502,t504,t52,t528,t53,t533,t54,t543,t544,t551,t552,t558,t56,t58,t59,t6,t60,t608,t62,t64,t67,t69,t7,t70,t72,t77,t78,t8,t80,t81,t82,t86,t91,t92,t95});
end
function fk_ub = ft_1(ct)
[L4,L5,L6,t10,t101,t102,t104,t11,t110,t113,t114,t115,t118,t12,t122,t129,t13,t130,t131,t133,t134,t135,t136,t137,t139,t14,t141,t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152,t153,t154,t155,t156,t16,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t195,t196,t197,t198,t199,t2,t201,t202,t205,t206,t209,t210,t211,t213,t214,t215,t216,t217,t218,t219,t221,t222,t223,t225,t226,t228,t229,t230,t231,t232,t233,t235,t238,t239,t240,t242,t243,t244,t245,t248,t250,t251,t252,t253,t254,t255,t257,t258,t260,t262,t264,t265,t266,t267,t268,t269,t270,t271,t272,t273,t274,t275,t276,t277,t278,t279,t280,t281,t282,t283,t284,t285,t286,t287,t288,t289,t290,t291,t292,t293,t294,t295,t296,t297,t298,t299,t300,t301,t302,t303,t304,t305,t306,t308,t309,t31,t310,t311,t312,t315,t316,t317,t319,t320,t322,t323,t324,t325,t328,t329,t330,t331,t334,t338,t339,t34,t340,t341,t343,t344,t345,t347,t348,t349,t35,t350,t351,t352,t353,t354,t355,t356,t357,t358,t359,t36,t360,t361,t362,t363,t364,t365,t366,t367,t369,t37,t370,t372,t373,t375,t376,t377,t378,t38,t380,t381,t382,t383,t4,t401,t402,t403,t404,t405,t406,t407,t409,t41,t410,t413,t414,t415,t416,t42,t420,t421,t422,t423,t428,t429,t43,t430,t434,t435,t436,t44,t441,t442,t447,t448,t45,t46,t469,t47,t470,t471,t472,t474,t48,t481,t483,t488,t49,t490,t495,t497,t5,t50,t502,t504,t52,t528,t53,t533,t54,t543,t544,t551,t552,t558,t56,t58,t59,t6,t60,t608,t62,t64,t67,t69,t7,t70,t72,t77,t78,t8,t80,t81,t82,t86,t91,t92,t95] = ct{:};
t611 = t150+t232+t255+t347+t410+t471;
t384 = -t381;
t385 = -t382;
t386 = -t383;
t387 = t381.*6.123233995736766e-17;
t389 = t383.*6.123233995736766e-17;
t391 = t381.*3.749399456654644e-33;
t393 = t383.*3.749399456654644e-33;
t394 = t381.*2.295845021658468e-49;
t396 = t383.*2.295845021658468e-49;
t397 = t381.*1.405799628556214e-65;
t399 = t383.*1.405799628556214e-65;
t417 = t47.*t414;
t418 = t52.*t414;
t424 = t416.*6.123233995736766e-17;
t431 = t416.*3.749399456654644e-33;
t437 = t416.*2.295845021658468e-49;
t443 = t416.*1.405799628556214e-65;
t449 = t34.*t447;
t450 = t41.*t447;
t451 = t47.*t448;
t452 = t52.*t448;
t473 = t139+t209+t265+t349;
t478 = t47.*t474;
t479 = t52.*t474;
t482 = -t481;
t484 = -t483;
t489 = -t488;
t496 = -t495;
t503 = -t502;
t530 = t45.*t528;
t531 = t50.*t528;
t591 = t45.*(-t144+t216+t229+t330+t420+t481);
t592 = t50.*(-t144+t216+t229+t330+t420+t481);
t609 = L4.*t608;
t612 = L4.*t611;
t388 = -t387;
t390 = -t389;
t392 = -t391;
t395 = -t394;
t398 = -t397;
t400 = -t399;
t419 = -t418;
t425 = -t424;
t426 = t418.*6.123233995736766e-17;
t432 = -t431;
t433 = t418.*3.749399456654644e-33;
t438 = -t437;
t439 = t418.*2.295845021658468e-49;
t444 = -t443;
t445 = t418.*1.405799628556214e-65;
t453 = -t451;
t454 = t449.*6.123233995736766e-17;
t456 = t451.*6.123233995736766e-17;
t458 = t449.*3.749399456654644e-33;
t460 = t451.*3.749399456654644e-33;
t461 = t449.*2.295845021658468e-49;
t463 = t451.*2.295845021658468e-49;
t464 = t449.*1.405799628556214e-65;
t465 = t451.*1.405799628556214e-65;
t475 = t34.*t473;
t476 = t41.*t473;
t480 = -t478;
t486 = t478.*6.123233995736766e-17;
t493 = t478.*3.749399456654644e-33;
t500 = t478.*2.295845021658468e-49;
t506 = t478.*1.405799628556214e-65;
t507 = t380+t450;
t508 = t385+t452;
t513 = -t37.*(t382-t452);
t514 = -t44.*(t382-t452);
t515 = t37.*(t382-t452);
t516 = t44.*(t382-t452);
t535 = t417+t479;
t541 = t531.*6.123233995736766e-17;
t549 = t531.*3.749399456654644e-33;
t556 = t531.*2.295845021658468e-49;
t575 = t135+t182+t185+t278+t384+t449;
t577 = t135+t182+t186+t282+t383+t451;
t593 = -t591;
t594 = -t37.*(t144+t215+t233+t334+t423+t484);
t596 = t37.*(t144+t215+t233+t334+t423+t484);
t597 = t591.*6.123233995736766e-17;
t610 = -t609;
t613 = t591.*3.749399456654644e-33;
t620 = t591.*2.295845021658468e-49;
t648 = t401+t612;
t654 = t530+t592;
t661 = -t12.*(t37.*(t409-t472)+t44.*(t144+t215+t233+t334+t423+t484));
t662 = t6.*(t37.*(t409-t472)+t44.*(t144+t215+t233+t334+t423+t484));
t665 = t12.*(t37.*(t409-t472)+t44.*(t144+t215+t233+t334+t423+t484)).*(-6.123233995736766e-17);
t672 = t12.*(t37.*(t409-t472)+t44.*(t144+t215+t233+t334+t423+t484)).*(-3.749399456654644e-33);
t427 = -t426;
t446 = -t445;
t455 = -t454;
t457 = -t456;
t466 = -t465;
t477 = -t475;
t485 = t475.*6.123233995736766e-17;
t487 = -t486;
t491 = t475.*3.749399456654644e-33;
t494 = -t493;
t498 = t475.*2.295845021658468e-49;
t501 = -t500;
t505 = t475.*1.405799628556214e-65;
t509 = t45.*t507;
t510 = t50.*t507;
t519 = t516.*(-6.123233995736766e-17);
t520 = t516.*6.123233995736766e-17;
t523 = t516.*(-3.749399456654644e-33);
t524 = t516.*3.749399456654644e-33;
t527 = t516.*(-2.295845021658468e-49);
t534 = t415+t476;
t539 = t37.*t535;
t540 = t44.*t535;
t542 = -t541;
t550 = -t549;
t576 = L4.*t575;
t578 = L4.*t577;
t598 = -t597;
t599 = t596.*(-6.123233995736766e-17);
t600 = t596.*6.123233995736766e-17;
t614 = -t613;
t615 = t596.*(-3.749399456654644e-33);
t621 = -t620;
t622 = t596.*(-2.295845021658468e-49);
t625 = L4.*(t151+t239-t265+t348-t416+t475);
t626 = t151+t239+t267+t350+t418+t480;
t647 = t402+t610;
t655 = t4.*t654;
t656 = t10.*t654;
t492 = -t491;
t512 = -t510;
t517 = t510.*6.123233995736766e-17;
t521 = t510.*3.749399456654644e-33;
t525 = t510.*2.295845021658468e-49;
t536 = t45.*t534;
t537 = t50.*t534;
t547 = t540.*6.123233995736766e-17;
t555 = t540.*3.749399456654644e-33;
t561 = t540.*2.295845021658468e-49;
t563 = t129+t170+t176+t264+t387+t455;
t564 = t129+t170+t174+t266+t390+t457;
t579 = -t578;
t587 = t376+t576;
t601 = t146+t217+t240+t341+t425+t485;
t602 = t145+t218+t242+t340+t427+t486;
t603 = -t45.*(t145+t218-t240+t338+t424-t485);
t604 = -t50.*(t145+t218-t240+t338+t424-t485);
t616 = t45.*(t145+t218-t240+t338+t424-t485).*(-6.123233995736766e-17);
t617 = t45.*(t145+t218-t240+t338+t424-t485).*6.123233995736766e-17;
t623 = t152+t238+t265+t349+t416+t477;
t627 = L4.*t626;
t629 = t45.*(t145+t218-t240+t338+t424-t485).*(-3.749399456654644e-33);
t633 = t45.*(t145+t218-t240+t338+t424-t485).*(-2.295845021658468e-49);
t652 = t406+t625;
t657 = -t655;
t658 = -t656;
t663 = t656.*6.123233995736766e-17;
t666 = t656.*3.749399456654644e-33;
t697 = t541+t598+t608;
t701 = t543+t600+t611;
t518 = -t517;
t522 = -t521;
t526 = -t525;
t538 = -t537;
t545 = t537.*6.123233995736766e-17;
t548 = -t547;
t553 = t537.*3.749399456654644e-33;
t559 = t537.*2.295845021658468e-49;
t562 = -t561;
t565 = t45.*t563;
t566 = t50.*t563;
t568 = t37.*t564;
t569 = t44.*t564;
t588 = t377+t579;
t605 = t37.*t602;
t606 = t44.*t602;
t628 = -t627;
t664 = -t663;
t667 = -t666;
t668 = t536+t604;
t698 = t4.*t697;
t699 = t10.*t697;
t702 = t6.*t701;
t703 = t12.*t701;
t546 = -t545;
t554 = -t553;
t560 = -t559;
t567 = -t565;
t570 = -t568;
t571 = t565.*6.123233995736766e-17;
t573 = t568.*6.123233995736766e-17;
t580 = t565.*3.749399456654644e-33;
t582 = t568.*3.749399456654644e-33;
t584 = t565.*2.295845021658468e-49;
t586 = t568.*2.295845021658468e-49;
t607 = -t605;
t618 = t605.*6.123233995736766e-17;
t631 = t605.*3.749399456654644e-33;
t634 = t605.*2.295845021658468e-49;
t638 = t10.*(t509-t566);
t639 = t515+t569;
t653 = t405+t628;
t669 = t4.*t668;
t670 = t10.*t668;
t673 = t539+t606;
t700 = -t698;
t704 = t698.*6.123233995736766e-17;
t706 = t702.*6.123233995736766e-17;
t707 = t698.*3.749399456654644e-33;
t709 = t702.*3.749399456654644e-33;
t711 = -t4.*(t151+t239-t265+t348-t416+t475+t545+t617);
t712 = -t10.*(t151+t239-t265+t348-t416+t475+t545+t617);
t713 = t4.*(t151+t239-t265+t348-t416+t475+t545+t617);
t737 = t657+t699;
t738 = -t35.*(t655-t699);
t739 = -t42.*(t655-t699);
t740 = t42.*(t655-t699);
t741 = t662+t703;
t788 = t165+t269+t286+t357+t428+t488+t542+t597+t656+t698;
t791 = t165+t269+t287+t359+t430+t490+t544+t599+t661+t702;
t572 = -t571;
t574 = -t573;
t581 = -t580;
t583 = -t582;
t585 = -t584;
t619 = -t618;
t632 = -t631;
t640 = t6.*t639;
t641 = t12.*t639;
t643 = t638.*(-6.123233995736766e-17);
t644 = t638.*6.123233995736766e-17;
t650 = t638.*3.749399456654644e-33;
t671 = -t670;
t674 = t6.*t673;
t675 = t12.*t673;
t677 = t670.*6.123233995736766e-17;
t681 = t670.*3.749399456654644e-33;
t683 = t517+t571+t575;
t687 = t519+t573+t577;
t705 = -t704;
t708 = -t707;
t710 = t546+t616+t623;
t714 = t548+t618+t626;
t717 = t713.*(-6.123233995736766e-17);
t718 = t713.*6.123233995736766e-17;
t722 = t713.*3.749399456654644e-33;
t742 = t48.*t741;
t743 = t53.*t741;
t745 = t740.*(-6.123233995736766e-17);
t748 = t669+t712;
t777 = t155+t253+t271+t352+t421+t482+t531+t593+t663+t704;
t781 = t155+t253+t274+t353+t423+t484+t533+t596+t665+t706;
et53 = t137.*2.295845021658468e-49+t251+t272-t353+t422+t483+t594-t706+t44.*(t409-t472);
et54 = t12.*(t37.*(t409-t472)+t44.*(t144+t215+t233+t334+t423+t484)).*6.123233995736766e-17;
t784 = -t48.*(et53+et54);
et55 = t137.*2.295845021658468e-49+t251+t272-t353+t422+t483+t594-t706+t44.*(t409-t472);
et56 = t12.*(t37.*(t409-t472)+t44.*(t144+t215+t233+t334+t423+t484)).*6.123233995736766e-17;
t785 = -t53.*(et55+et56);
et57 = t137.*2.295845021658468e-49+t251+t272-t353+t422+t483+t594-t706+t44.*(t409-t472);
et58 = t12.*(t37.*(t409-t472)+t44.*(t144+t215+t233+t334+t423+t484)).*6.123233995736766e-17;
t786 = t48.*(et57+et58).*(-6.123233995736766e-17);
et59 = t137.*2.295845021658468e-49+t251+t272-t353+t422+t483+t594-t706+t44.*(t409-t472);
et60 = t12.*(t37.*(t409-t472)+t44.*(t144+t215+t233+t334+t423+t484)).*6.123233995736766e-17;
t787 = t48.*(et59+et60).*6.123233995736766e-17;
t789 = L5.*t788;
t792 = L5.*t791;
t806 = L5.*(t166-t273-t291+t360-t431+t491+t546+t616+t670+t713);
t642 = -t641;
t645 = t641.*6.123233995736766e-17;
t651 = t641.*3.749399456654644e-33;
t676 = -t675;
t678 = -t677;
t679 = t675.*6.123233995736766e-17;
t682 = t675.*3.749399456654644e-33;
t684 = t4.*t683;
t685 = t10.*t683;
t688 = t6.*t687;
t689 = t12.*t687;
t715 = t6.*t714;
t716 = t12.*t714;
t744 = -t743;
t746 = t743.*6.123233995736766e-17;
t749 = t35.*t748;
t750 = t42.*t748;
t778 = t35.*t777;
t779 = t42.*t777;
t790 = -t789;
t794 = -t42.*(-t156+t257+t277-t354+t424-t485+t537+t677+t718+t45.*(t145+t218-t240+t338+t424-t485));
t796 = t35.*(-t156+t257+t277-t354+t424-t485+t537+t677+t718+t45.*(t145+t218-t240+t338+t424-t485));
t804 = t167+t273+t291+t361+t431+t492+t545+t617+t671+t711;
t811 = t648+t792;
t812 = t652+t806;
t830 = t742+t785;
t646 = -t645;
t680 = -t679;
t690 = t684.*6.123233995736766e-17;
t692 = t688.*6.123233995736766e-17;
t694 = t684.*3.749399456654644e-33;
t695 = t688.*3.749399456654644e-33;
t719 = t715.*6.123233995736766e-17;
t723 = t715.*3.749399456654644e-33;
t726 = t35.*(t685-t4.*(t509-t566));
t727 = t42.*(t685-t4.*(t509-t566));
t729 = t640+t689;
t747 = -t746;
t752 = t674+t716;
t757 = t750.*6.123233995736766e-17;
t771 = L5.*(t147+t197+t201+t308+t392+t458+t518+t572+t638+t684);
t772 = t147+t197+t202+t310+t393+t460+t520+t574+t642+t688;
t780 = -t778;
t782 = t778.*6.123233995736766e-17;
t793 = t156+t258+t280+t354+t425+t485+t538+t603+t678+t717;
t798 = t796.*(-6.123233995736766e-17);
t799 = t796.*6.123233995736766e-17;
t807 = t166+t275+t295+t362+t433+t494+t547+t619+t676+t715;
t810 = t647+t790;
t826 = t738+t779;
t827 = t36.*(t779-t35.*(t655-t699));
t828 = t43.*(t779-t35.*(t655-t699));
t831 = t49.*t830;
t832 = t54.*t830;
t835 = t749+t794;
t859 = t746+t787+t791;
t691 = -t690;
t693 = -t692;
t696 = -t695;
t720 = -t719;
t724 = -t723;
t728 = -t727;
t730 = t48.*t729;
t731 = t53.*t729;
t734 = t727.*6.123233995736766e-17;
t753 = t48.*t752;
t754 = t53.*t752;
t759 = t142+t193+t194+t290+t387+t455+t510+t565+t644+t690;
t773 = L5.*t772;
t775 = t587+t771;
t783 = -t782;
t808 = L5.*t807;
t829 = t828.*6.123233995736766e-17;
t833 = -t832;
t834 = t832.*6.123233995736766e-17;
t836 = t36.*t835;
t837 = t43.*t835;
t860 = t49.*t859;
t861 = t54.*t859;
t865 = t757+t799+t804;
t732 = -t730;
t733 = -t731;
t735 = -t734;
t736 = t731.*6.123233995736766e-17;
t755 = -t753;
t756 = -t754;
t758 = t754.*6.123233995736766e-17;
t760 = t35.*t759;
t761 = t42.*t759;
t763 = t141+t192+t196+t296+t389+t456+t516+t570+t645+t693;
t774 = -t773;
t797 = t156+t258+t281+t355+t426+t487+t540+t607+t679+t720;
t809 = -t808;
t839 = t837.*6.123233995736766e-17;
t854 = t745+t783+t788;
t862 = -t860;
t863 = -t861;
t864 = t860.*6.123233995736766e-17;
t866 = t36.*t865;
t867 = t43.*t865;
t764 = t48.*t763;
t765 = t53.*t763;
t767 = t760.*6.123233995736766e-17;
t776 = t588+t774;
t800 = t53.*t797;
t801 = t48.*t797;
t813 = t653+t809;
t814 = t726+t761;
t855 = t36.*t854;
t856 = t43.*t854;
t868 = t866.*6.123233995736766e-17;
t876 = t831+t863;
t884 = t179+t283+t301+t366+t436+t497+t552+t615+t665+t706+t744+t784+t834+t864;
t886 = t189+t297+t312+t372+t442+t504+t558+t622+t672+t709+t747+t786+t833+t862;
t889 = t190+t303+t317+t373+t444+t505+t560+t633+t681+t722+t757+t799+t837+t866;
t766 = -t764;
t768 = t764.*6.123233995736766e-17;
t802 = -t801;
t803 = t801.*6.123233995736766e-17;
t815 = t36.*t814;
t816 = t43.*t814;
t820 = t732+t765;
t821 = -t49.*(t730-t765);
t823 = t54.*(t730-t765);
t840 = t755+t800;
t841 = -t49.*(t753-t800);
t843 = t54.*(t753-t800);
t847 = -t36.*(t147+t197+t201+t308+t392+t458+t518+t572+t638+t684+t734-t767);
t848 = -t43.*(t147+t197+t201+t308+t392+t458+t518+t572+t638+t684+t734-t767);
t857 = -t856;
t858 = t855.*6.123233995736766e-17;
t885 = t188+t298+t311+t370+t441+t503+t556+t621+t667+t708+t745+t783+t828+t855;
t817 = -t816;
t818 = t816.*6.123233995736766e-17;
t825 = t823.*6.123233995736766e-17;
t844 = t843.*(-6.123233995736766e-17);
t845 = t843.*6.123233995736766e-17;
t850 = t736+t768+t772;
t869 = t758+t803+t807;
t873 = t815+t848;
t875 = t827+t857;
t883 = t179+t283+t300+t364+t434+t495+t550+t613+t663+t704+t740+t778+t829+t858;
t851 = t49.*t850;
t852 = t54.*t850;
t870 = t49.*t869;
t871 = t54.*t869;
t881 = t162+t221+t226+t322+t398+t464+t526+t585+t650+t694+t735+t767+t817+t847;
t853 = t851.*6.123233995736766e-17;
t872 = t870.*6.123233995736766e-17;
t874 = t821+t852;
t878 = t841+t871;
t882 = t164+t222+t228+t323+t400+t466+t527+t586+t651+t696+t736+t768+t823+t851;
t890 = t191+t302+t320+t375+t446+t506+t562+t634+t682+t724+t758+t803+t843+t870;
t880 = t153+t210+t214+t319+t396+t463+t524+t583+t646+t692+t733+t766+t825+t853;
t888 = t180+t289+t306+t369+t439+t501+t555+t632+t680+t719+t756+t802+t845+t872;
et61 = t137.*2.295845021658468e-49+t251+t272-t353+t422+t483+t594-t706+t44.*(t409-t472);
et62 = t12.*(t37.*(t409-t472)+t44.*(t144+t215+t233+t334+t423+t484)).*6.123233995736766e-17;
et63 = t178+t284-t301+t365-t436-t497+t551+t596.*3.749399456654644e-33-t706+t743-t834-t864+t12.*(t37.*(t409-t472)+t44.*(t144+t215+t233+t334+t423+t484)).*6.123233995736766e-17;
et64 = t48.*(et61+et62);
fk_ub = reshape([t2,0.0,t16,0.0,0.0,1.0,0.0,0.0,t8,0.0,t2,0.0,0.0,0.0,0.0,1.0,t56,1.0,t58,0.0,t8-t56,3.749399456654644e-33,t59,0.0,t59,-6.123233995736766e-17,t60,0.0,t15,0.0,t14,1.0,t95,t78,t69-t8.*t31.*6.123233995736766e-17,0.0,t114,t62+t92-6.123233995736766e-17,t113,0.0,t56+t64+t72+t82,t38+t77-3.749399456654644e-33,t58+t67+t70+t81,0.0,t15,0.0,t14,1.0,t133,t122,t131,0.0,t206,t38+t77-t102-t115-3.749399456654644e-33,t199,0.0,t231,t168,t225,0.0,t15,0.0,t14,1.0,t328,t248,t324,0.0,t331,t250,t325,0.0,t231,t168,t225,0.0,t15+t235,t169,t14+t223,1.0,t328,t248,t324,0.0,t231,t168,t225,0.0,t139+t209,t245,t134+t205,0.0,t406,t376,t402,1.0,t413,t378,t403,0.0,t473,t447,t134+t205+t254+t345,0.0,t145+t218-t240+t338,t129+t170+t176+t264,t143+t216+t229+t330,0.0,t406,t376,t402,1.0,t534,t507,t407-t470,0.0,t601,t130+t171+t172+t260+t388+t454,t144+t215-t229-t330+t421+t482,0.0,t151+t239-t265+t348-t416+t475,t575,t149-t232+t252+t344-t469+t41.*(t230-t329),0.0,t406,t376,t402,1.0,t668,t509-t566,-t654,0.0,t710,t136+t183+t184+t276+t381-t449+t518+t572,t697,0.0,t156+t258+t280+t354+t425+t485+t538+t603,t141+t192+t195+t292+t388+t454+t512+t567,t155+t253+t271+t352+t421+t482+t531+t593,0.0,t652,t587,t647,1.0,t748,-t685+t4.*(t509-t566),t737,0.0,t793,t141+t192+t195+t292+t388+t454+t512+t567+t643+t691,t777,0.0,t166+t275+t294+t360+t432+t491+t546+t616+t670+t713,t147+t197+t201+t308+t392+t458+t518+t572+t638+t684,t163+t270+t285+t356+t429+t489+t541+t598+t658+t700,0.0,t652,t587,t647,1.0,t835,-t814,t826,0.0,t166+t275+t294+t360+t432+t491+t546+t616+t670+t713-t757+t798,t147+t197+t201+t308+t392+t458+t518+t572+t638+t684+t734-t767,t163+t270+t285+t356+t429+t489+t541+t598+t658+t700+t740.*6.123233995736766e-17+t782,0.0,t180+t289+t305+t367+t438+t498+t554+t629+t677+t718+t750+t796,t153+t210+t213+t315+t395+t461+t522+t581+t644+t690+t728+t760,t178+t284+t299+t363+t435+t496+t549+t614+t664+t705+t739+t780,0.0,t812,t775,t810,1.0,t836-t867,-t815+t43.*(t147+t197+t201+t308+t392+t458+t518+t572+t638+t684+t734-t767),t875,0.0,t180+t289+t305+t367+t438+t498+t554+t629+t677+t718+t750+t796-t839-t868,t153+t210+t213+t315+t395+t461+t522+t581+t644+t690+t728+t760+t818+t36.*(t147+t197+t201+t308+t392+t458+t518+t572+t638+t684+t734-t767).*6.123233995736766e-17,t178+t284+t299+t363+t435+t496+t549+t614+t664+t705+t739+t780-t829-t858,0.0,t889,t881,t885,0.0,t812,t775,t810,1.0,t11.*(t180-t288-t304+t367-t437+t498+t554+t629+t677+t718+t750+t796-t839-t868)+t5.*(t836-t867),-t5.*t873+t11.*(t153+t210+t213+t315+t395+t461+t522+t581+t644+t690+t728+t760+t818+t36.*(t147+t197+t201+t308+t392+t458+t518+t572+t638+t684+t734-t767).*6.123233995736766e-17),t5.*t875-t11.*t883,0.0,t5.*(t180-t288-t304+t367-t437+t498+t554+t629+t677+t718+t750+t796-t839-t868)-t11.*(t836-t867),t11.*t873+t5.*(t153+t210+t213+t315+t395+t461+t522+t581+t644+t690+t728+t760+t818+t36.*(t147+t197+t201+t308+t392+t458+t518+t572+t638+t684+t734-t767).*6.123233995736766e-17),-t11.*t875-t5.*t883,0.0,t889,t881,t885,0.0,t812+L6.*t11.*(t180-t288-t304+t367-t437+t498+t554+t629+t677+t718+t750+t796-t839-t868)+L6.*t5.*(t836-t867),t775-L6.*t5.*t873+L6.*t11.*(t153+t210+t213+t315+t395+t461+t522+t581+t644+t690+t728+t760+t818+t36.*(t147+t197+t201+t308+t392+t458+t518+t572+t638+t684+t734-t767).*6.123233995736766e-17),t810+L6.*t5.*t875-L6.*t11.*t883,1.0,t328,t248,t324,0.0,-t81+t86+t91-t101-t104+t118,t62-t80+t92-t110+2.295845021658468e-49,t219,0.0,t331,t250,t325,0.0,t405,t377,t401,1.0,-t243+t339,-t175+t46.*t248,t404,0.0,t139+t209+t268+t351,t448,t134+t205+t255+t347,0.0,t146+t217+t244+t343,t130+t171+t177+t262,t144+t215+t233+t334,0.0,t405,t377,t401,1.0,-t535,t508,t409-t472,0.0,t146+t217+t244+t343+t426+t487,t130+t171+t177+t262+t389+t456,t144+t215+t233+t334+t423+t484,0.0,t152+t238+t268+t351+t419+t478,t136+t183+t187+t279+t386+t453,t611,0.0,t405,t377,t401,1.0,-t673,t513-t569,t37.*(t409-t472)+t44.*(t144+t215+t233+t334+t423+t484),0.0,t152+t238+t268+t351+t419+t478+t547+t619,t136+t183+t187+t279+t386+t453+t520+t574,t701,0.0,-t156+t257-t281-t355+t427+t486-t540+t605,t142+t193-t196+t293+t390+t457+t514+t568,t137.*2.295845021658468e-49+t251+t272-t353+t422+t483+t594+t44.*(t409-t472),0.0,t653,t588,t648,1.0,-t752,-t729,t741,0.0,t797,t763,t781,0.0,t167+t273-t295-t362-t433+t493+t548+t618+t675-t715,t148+t198-t202+t309-t393-t460+t519+t573+t641-t688,t791,0.0,t653,t588,t648,1.0,t840,t820,t830,0.0,t869,t850,t163+t270-t287+t358-t430-t490+t543+t600-t702+t747+t786+t12.*(t37.*(t409-t472)+t44.*(t144+t215+t233+t334+t423+t484)),0.0,t181+t288-t306-t369-t439+t500-t555+t631+t679+t720+t754+t801,t154+t211-t214+t316-t396-t463+t523+t582+t645+t693+t731+t764,t179+t283+t301+t366+t436+t497+t552+t615+t665+t706+t744+t784,0.0,t813,t776,t811,1.0,t878,t874,t876,0.0,t888,t880,et63+et64,0.0,t890,t882,t886,0.0,t813,t776,t811,1.0,t7.*(t871-t49.*(t753-t800))-t13.*(t181+t288-t306-t369-t439+t500-t555+t631+t679+t720+t754+t801+t844-t872),t7.*(t852-t49.*(t730-t765))+t13.*t880,t7.*t876-t13.*t884,0.0,-t13.*(t871-t49.*(t753-t800))-t7.*(t181+t288-t306-t369-t439+t500-t555+t631+t679+t720+t754+t801+t844-t872),-t13.*(t852-t49.*(t730-t765))+t7.*t880,-t13.*t876-t7.*t884,0.0,t890,t882,t886,0.0,t813+L6.*t7.*(t871-t49.*(t753-t800))-L6.*t13.*(t181+t288-t306-t369-t439+t500-t555+t631+t679+t720+t754+t801+t844-t872),t776+L6.*t7.*(t852-t49.*(t730-t765))+L6.*t13.*t880,t811+L6.*t7.*t876-L6.*t13.*t884,1.0],[4,4,21]);
end

function dJ_ll_dq = dJ_ll_dq_computable(in1,in2,in3)
%dJ_ll_dq_computable
%    dJ_ll_dq = dJ_ll_dq_computable(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    22-Mar-2024 14:10:13

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
t2 = cos(ql5);
t3 = cos(ql6);
t4 = sin(ql5);
t5 = sin(ql6);
t6 = Ll5+Ll6;
t7 = -xg5;
t8 = -xg8;
t9 = -xg11;
t10 = pi./2.0;
t11 = -t10;
t12 = ql3+t10;
t13 = ql4+t10;
t14 = ql1+t11;
t15 = ql2+t11;
t16 = ql7+t11;
t17 = cos(t12);
t18 = cos(t13);
t19 = sin(t12);
t20 = sin(t13);
t21 = cos(t14);
t22 = cos(t15);
t23 = cos(t16);
t24 = sin(t14);
t25 = sin(t15);
t26 = sin(t16);
t27 = Ll2.*t25;
t28 = -t25;
t29 = t17.*t22;
t30 = t19.*t22;
t31 = t17.*t25;
t32 = t19.*t25;
t33 = t21.*t22;
t34 = t21.*t25;
t35 = t22.*t24;
t36 = t24.*t25;
t40 = t21.*6.123233995736766e-17;
t41 = t22.*6.123233995736766e-17;
t42 = t24.*6.123233995736766e-17;
t43 = t25.*6.123233995736766e-17;
t76 = t21.*3.749399456654644e-33;
t77 = t22.*3.749399456654644e-33;
t78 = t22+3.749399456654644e-33;
t79 = t24.*3.749399456654644e-33;
t80 = t25.*3.749399456654644e-33;
t122 = t21.*2.295845021658468e-49;
t123 = t24.*2.295845021658468e-49;
t37 = t19.*t28;
t38 = t21.*t28;
t39 = t24.*t28;
t44 = -t42;
t45 = t30.*6.123233995736766e-17;
t46 = t31.*6.123233995736766e-17;
t47 = t32.*6.123233995736766e-17;
t48 = t33.*6.123233995736766e-17;
t49 = t34.*6.123233995736766e-17;
t50 = t35.*6.123233995736766e-17;
t51 = t36.*6.123233995736766e-17;
t83 = -t79;
t84 = Ll2.*t78;
t85 = t41-6.123233995736766e-17;
t86 = t31.*3.749399456654644e-33;
t87 = t33.*3.749399456654644e-33;
t88 = t34.*3.749399456654644e-33;
t89 = t35.*3.749399456654644e-33;
t90 = t36.*3.749399456654644e-33;
t124 = -t123;
t125 = t33.*2.295845021658468e-49;
t126 = t34.*2.295845021658468e-49;
t127 = t35.*2.295845021658468e-49;
t128 = t36.*2.295845021658468e-49;
t52 = -t47;
t53 = -t48;
t54 = -t50;
t55 = -t51;
t57 = t34+t50;
t58 = t35+t49;
t60 = t39+t48;
t91 = -t87;
t93 = t17.*t85;
t94 = t19.*t85;
t129 = -t125;
t133 = t50+t88;
t135 = t24+t49+t89;
t141 = t28+t45+t86;
t56 = t29+t52;
t59 = t33+t55;
t61 = Ll2.*t58;
t65 = t17.*t57;
t66 = t17.*t58;
t67 = t19.*t57;
t68 = t19.*t58;
t72 = t17.*t60;
t73 = t19.*t60;
t95 = -t93;
t96 = t31+t94;
t106 = t36+t40+t53;
t109 = t44+t57;
t130 = t93.*6.123233995736766e-17;
t131 = t94.*6.123233995736766e-17;
t136 = t17.*t133;
t137 = t19.*t133;
t138 = t21+t55+t87;
t140 = -t19.*(t53+t90);
t143 = t17.*t135;
t144 = t19.*t135;
t147 = t18.*t141;
t148 = t20.*t141;
t158 = t17.*(t53+t90).*6.123233995736766e-17;
t62 = Ll2.*t59;
t63 = t18.*t56;
t64 = t20.*t56;
t70 = t17.*t59;
t71 = t19.*t59;
t97 = t32+t95;
t98 = t66.*6.123233995736766e-17;
t99 = t67.*6.123233995736766e-17;
t100 = t68.*6.123233995736766e-17;
t101 = t18.*t96;
t102 = t20.*t96;
t110 = t73.*6.123233995736766e-17;
t115 = Ll2.*t106;
t116 = Ll2.*t109;
t132 = -t130;
t142 = -t137;
t145 = t17.*t138;
t146 = t19.*t138;
t149 = Ll3.*t148;
t150 = -t144;
t152 = -t148;
t154 = t46+t131;
t155 = t136.*6.123233995736766e-17;
t163 = t143.*6.123233995736766e-17;
t164 = t144.*6.123233995736766e-17;
t169 = t65+t140;
t69 = Ll3.*t63;
t103 = t70.*6.123233995736766e-17;
t104 = t71.*6.123233995736766e-17;
t105 = -t100;
t108 = Ll3.*t101;
t111 = Ll3.*t102;
t112 = t18.*t97;
t113 = t20.*t97;
t120 = t115+xg3;
t151 = -t145;
t153 = -t149;
t159 = t18.*t154;
t160 = t20.*t154;
t165 = t145.*6.123233995736766e-17;
t166 = t146.*6.123233995736766e-17;
t167 = -t164;
t170 = t72+t142;
t171 = t64+t147;
t172 = t66+t146;
t173 = t71+t143;
t174 = t18.*t169;
t175 = t20.*t169;
t176 = t63+t152;
t179 = t70+t150;
t207 = t47+t78+t132;
t238 = t59+t99+t158;
t243 = -t20.*(-t58+t110+t155);
t107 = -t103;
t117 = Ll3.*t112;
t118 = -t108;
t161 = Ll3.*t160;
t162 = -t159;
t168 = -t165;
t177 = Ll3.*t174;
t178 = t68+t151;
t180 = t2.*t171;
t181 = t4.*t171;
t182 = t18.*t170;
t183 = t20.*t170;
t185 = t2.*t176;
t187 = t4.*t176;
t190 = t18.*t172;
t191 = t18.*t173;
t192 = t20.*t172;
t193 = t20.*t173;
t198 = t18.*t179;
t200 = t20.*t179;
t208 = t18.*t207;
t209 = t20.*t207;
t217 = t112+t160;
t221 = t98+t166;
t225 = t103+t167;
t240 = t18.*t238;
t241 = t20.*t238;
t246 = Ll3.*t243;
t257 = t105+t106+t165;
t258 = t38+t42+t54+t104+t163;
t184 = -t177;
t186 = Ll3.*t182;
t194 = Ll3.*t190;
t195 = Ll3.*t191;
t196 = Ll3.*t192;
t197 = t18.*t178;
t199 = t20.*t178;
t202 = Ll3.*t198;
t204 = Ll3.*t200;
t210 = t180.*6.123233995736766e-17;
t211 = Ll3.*t208;
t212 = Ll3.*t209;
t214 = t187.*6.123233995736766e-17;
t215 = -t209;
t218 = t113+t162;
t219 = t2.*t217;
t220 = t4.*t217;
t228 = t18.*t221;
t229 = t20.*t221;
t231 = t18.*t225;
t232 = t20.*t225;
t244 = Ll3.*t241;
t245 = -t241;
t247 = t102+t208;
t264 = t18.*t257;
t265 = t20.*t257;
t266 = t18.*t258;
t267 = t20.*t258;
t275 = -Ll4.*t3.*(t181-t185);
t297 = t175+t240;
t301 = t182+t243;
t309 = -t4.*(t183+t18.*(-t58+t110+t155));
t310 = t2.*(t183+t18.*(-t58+t110+t155));
t334 = t45+t80+t86+t180+t187;
t201 = Ll3.*t197;
t213 = -t210;
t216 = -t214;
t222 = t2.*t218;
t223 = t4.*t218;
t230 = Ll3.*t229;
t233 = Ll3.*t232;
t234 = -t228;
t235 = -t231;
t236 = t220.*6.123233995736766e-17;
t248 = t101+t215;
t249 = t2.*t247;
t250 = t4.*t247;
t268 = Ll3.*t264;
t269 = Ll3.*t265;
t270 = Ll3.*t266;
t271 = Ll3.*t267;
t272 = -t264;
t276 = Ll1+t7+t84+t118+t212;
t277 = t191+t232;
t278 = t197+t229;
t298 = t174+t245;
t299 = t2.*t297;
t300 = t4.*t297;
t306 = t2.*t301;
t307 = t4.*t301;
t317 = t310.*(-6.123233995736766e-17);
t320 = t190+t265;
t323 = t200+t266;
t331 = -t2.*(t198-t267);
t333 = t4.*(t198-t267);
t335 = Ll5.*t334;
t345 = t2.*(t198-t267).*(-6.123233995736766e-17);
t227 = -t223;
t237 = t222.*6.123233995736766e-17;
t251 = t2.*t248;
t252 = -t249;
t253 = t4.*t248;
t259 = t249.*6.123233995736766e-17;
t260 = t250.*6.123233995736766e-17;
t273 = -t271;
t279 = t193+t235;
t280 = t2.*t277;
t281 = t4.*t277;
t282 = t199+t234;
t283 = t2.*t278;
t284 = t4.*t278;
t302 = t2.*t298;
t303 = t4.*t298;
t313 = t299.*6.123233995736766e-17;
t315 = t307.*6.123233995736766e-17;
t318 = t30+t43+t46+t213+t216;
t321 = t192+t272;
t322 = t4.*t320;
t324 = t2.*t320;
t329 = t2.*t323;
t330 = t4.*t323;
t336 = -t335;
t346 = t333.*(-6.123233995736766e-17);
t359 = t115+t194+t269+xg6;
t364 = -Ll5.*(-t154+t220+t222);
t389 = t306+t309;
t451 = t89+t110+t126+t155+t307+t310;
t255 = -t251;
t256 = -t253;
t261 = -t260;
t262 = t251.*6.123233995736766e-17;
t263 = t253.*6.123233995736766e-17;
t285 = t2.*t279;
t286 = t4.*t279;
t287 = t2.*t282;
t288 = t4.*t282;
t293 = t281.*6.123233995736766e-17;
t294 = t284.*6.123233995736766e-17;
t311 = t219+t227;
t314 = t303.*6.123233995736766e-17;
t316 = -t315;
t319 = Ll4.*t5.*t318;
t325 = t4.*t321;
t327 = t2.*t321;
t337 = t324.*6.123233995736766e-17;
t338 = t322.*6.123233995736766e-17;
t342 = t329.*6.123233995736766e-17;
t343 = t330.*6.123233995736766e-17;
t347 = t249+t253;
t361 = t96+t236+t237;
t388 = Ll4.*t3.*(t300-t302);
t390 = Ll4.*t3.*t389;
t399 = t330+t331;
t405 = t329+t333;
t449 = t91+t99+t128+t158+t299+t303;
t452 = Ll5.*t451;
t290 = -t286;
t292 = -t288;
t295 = t285.*6.123233995736766e-17;
t296 = t287.*6.123233995736766e-17;
t312 = Ll4.*t3.*t311;
t339 = -t337;
t340 = t327.*6.123233995736766e-17;
t341 = t325.*6.123233995736766e-17;
t344 = -t342;
t349 = t3.*t347;
t351 = -Ll5.*(t250+t255);
t356 = t3.*(t250+t255);
t357 = t6.*(t250+t255);
t363 = Ll4.*t5.*t361;
t366 = t5.*t347.*6.123233995736766e-17;
t368 = t261+t262;
t373 = t37+t41+t93+t259+t263+2.295845021658468e-49;
t382 = t52+t77+t130+t252+t256+1.405799628556214e-65;
t391 = t322+t327;
t398 = -t6.*(t324-t325);
t403 = Ll5.*t399;
t406 = Ll4.*t3.*(t324-t325);
t407 = t3.*t399;
t408 = t6.*t399;
t410 = Ll4.*t5.*t399;
t412 = t3.*t405;
t419 = t5.*(t324-t325).*6.123233995736766e-17;
t420 = t5.*t399.*6.123233995736766e-17;
t422 = t5.*t405.*6.123233995736766e-17;
t426 = t343+t345;
t427 = t107+t164+t281+t285;
t438 = -Ll5.*(-t221+t284+t287);
t446 = -Ll4.*t5.*(t53+t67+t90-t313-t314+t17.*(t53+t90));
t447 = t73+t133+t136+t316+t317;
t450 = Ll5.*t449;
t453 = -t452;
t473 = t88+t104+t124+t127+t163+t405;
t350 = Ll4.*t349;
t358 = Ll4.*t356;
t365 = -t363;
t369 = t5.*t368;
t371 = t3.*t368.*6.123233995736766e-17;
t374 = t280+t290;
t375 = t5.*t373;
t376 = t283+t292;
t377 = Ll4.*t3.*t373;
t381 = t3.*t373.*6.123233995736766e-17;
t383 = Ll5.*t382;
t384 = t6.*t382;
t393 = t3.*t391;
t409 = Ll4.*t407;
t411 = -t403;
t414 = Ll4.*t412;
t417 = t5.*t391.*6.123233995736766e-17;
t421 = -t420;
t423 = t339+t341;
t428 = t172+t294+t296;
t429 = t179+t293+t295;
t430 = Ll5.*t427;
t431 = t5.*t426;
t442 = t3.*t426.*6.123233995736766e-17;
t448 = Ll4.*t5.*t447;
t455 = t5.*(-t51-t68-t76+t87+t145+t338+t340);
t456 = t49+t83+t89+t173+t344+t346;
t457 = Ll4.*t3.*(-t51-t68-t76+t87+t145+t338+t340);
t466 = t3.*(-t51-t68-t76+t87+t145+t338+t340).*6.123233995736766e-17;
t469 = t90+t100+t122+t129+t168+t391;
t474 = Ll5.*t473;
t475 = t6.*t473;
t370 = Ll4.*t369;
t372 = -t371;
t378 = Ll4.*t375;
t379 = Ll4.*t3.*t374;
t380 = Ll4.*t3.*t376;
t385 = -t383;
t395 = Ll4.*t393;
t416 = -t414;
t424 = t5.*t423;
t432 = Ll4.*t431;
t433 = -t430;
t434 = Ll4.*t5.*t429;
t437 = Ll4.*t5.*t428;
t441 = t3.*t423.*6.123233995736766e-17;
t443 = t349+t369;
t458 = Ll4.*t455;
t460 = t5.*t456;
t461 = Ll4.*t3.*t456;
t464 = t356+t375;
t467 = t3.*t456.*6.123233995736766e-17;
t470 = Ll5.*t469;
t471 = t6.*t469;
t476 = -t475;
t483 = -Ll7.*t23.*(t412-t431);
t485 = -Ll7.*t26.*(-t381+t382+t5.*(t250+t255).*6.123233995736766e-17);
t488 = Ll7.*t23.*(t455-t3.*(t324-t325));
t495 = t399+t422+t442;
t497 = t419+t466+t469;
t425 = Ll4.*t424;
t436 = -t432;
t439 = -t434;
t440 = -t437;
t444 = Ll7.*t23.*t443;
t459 = -t458;
t462 = Ll4.*t460;
t465 = Ll7.*t23.*t464;
t468 = -t467;
t472 = -t470;
t477 = t250+t255+t366+t372;
t489 = -t488;
t491 = Ll7.*t23.*(t407-t460);
t492 = Ll1+t8+t84+t118+t212+t358+t378+t385;
t494 = Ll7.*t26.*(-t324+t325+t417+t441);
t496 = Ll7.*t26.*t495;
t498 = Ll7.*t26.*t497;
t463 = -t462;
t478 = Ll7.*t26.*t477;
t499 = -t498;
t500 = t421+t468+t473;
t502 = t115+t194+t269+t406+t459+t472+xg9;
t504 = Ll1+t9+t84+t118+t212+t358+t378+t384+t385+t465+t485;
t479 = -t478;
t501 = Ll7.*t26.*t500;
t505 = t115+t194+t269+t406+t459+t471+t472+t489+t499+xg12;
et1 = t62-t177+t244+t388+t446+t450-t6.*t449-Ll7.*t23.*(t5.*(t53+t67+t90-t313-t314+t17.*(t53+t90))-t3.*(t300-t302))-Ll7.*t26.*(-t449+t3.*(t53+t67+t90-t313-t314+t17.*(t53+t90)).*6.123233995736766e-17+t5.*(t300-t302).*6.123233995736766e-17);
et2 = (-t116-t202+t271+t409+t463+t474+t476+t491+t501+xg10).*-2.0;
et3 = t505.*(t61+t186+t246+t390+t448+t453+t6.*t451+Ll7.*t23.*(t3.*t389+t5.*t447)-Ll7.*t26.*(t451+t5.*t389.*6.123233995736766e-17-t3.*t447.*6.123233995736766e-17)).*2.0+t61.*t120.*2.0+t62.*(t116-xg1).*2.0-t276.*(t27+t69+t153).*2.0+t359.*(t61+t186+t246).*2.0+et1.*et2+(t62+t184+t244).*(t116+t202+t273-xg4).*2.0;
et4 = (t62+t184+t244+t388+t446+t450).*(-t116-t202+t271+t409+t463+t474+xg7).*-2.0-t27.*(Ll1+t84-xg2).*2.0-t492.*(t27+t69+t153+t275+t319+t336).*2.0+t502.*(t61+t186+t246+t390+t448+t453).*2.0;
et5 = t504.*(t27+t69+t153+t275+t319+t336+t6.*t334+Ll7.*t23.*(t5.*t318-t3.*(t181-t185))+Ll7.*t26.*(-t334+t3.*t318.*6.123233995736766e-17+t5.*(t181-t185).*6.123233995736766e-17)).*-2.0;
et6 = t276.*(t117+t161).*2.0-t359.*(t201+t230).*2.0-(t195+t233).*(t116+t202+t273-xg4).*2.0+(t195+t233+t379+t433+t439).*(-t116-t202+t271+t409+t463+t474+xg7).*2.0+t492.*(t117+t161+t312+t364+t365).*2.0-t502.*(t201+t230+t380+t438+t440).*2.0;
et7 = t504.*(t117+t161+t312+t364+t365+t6.*(-t154+t220+t222)+Ll7.*t23.*(t3.*t311-t5.*t361)-Ll7.*t26.*(-t154+t220+t222+t5.*t311.*6.123233995736766e-17+t3.*t361.*6.123233995736766e-17)).*2.0;
et8 = t505.*(t201+t230+t380+t438+t440+t6.*(-t221+t284+t287)+Ll7.*t23.*(t3.*t376-t5.*t428)-Ll7.*t26.*(-t221+t284+t287+t5.*t376.*6.123233995736766e-17+t3.*t428.*6.123233995736766e-17)).*-2.0;
et9 = (t195+t233+t379+t433+t439+t6.*t427+Ll7.*t23.*(t3.*t374-t5.*t429)-Ll7.*t26.*(t427+t5.*t374.*6.123233995736766e-17+t3.*t429.*6.123233995736766e-17)).*(-t116-t202+t271+t409+t463+t474+t476+t491+t501+xg10).*2.0;
et10 = t492.*(t377-Ll4.*t5.*(t250+t255)).*2.0-t505.*(t457-Ll7.*t26.*(t455.*6.123233995736766e-17-t3.*(t324-t325).*6.123233995736766e-17)+Ll7.*t23.*(t5.*(t324-t325)+t3.*(-t51-t68-t76+t87+t145+t338+t340))+Ll4.*t5.*(t324-t325)).*2.0-t502.*(t457+Ll4.*t5.*(t324-t325)).*2.0;
et11 = t504.*(-t377+Ll4.*t5.*(t250+t255)+Ll7.*t26.*(t356.*6.123233995736766e-17+t375.*6.123233995736766e-17)+Ll7.*t23.*(t5.*(t250+t255)-t3.*t373)).*-2.0;
et12 = (t410+t461+Ll7.*t23.*(t5.*t399+t3.*t456)+Ll7.*t26.*(t407.*6.123233995736766e-17-t460.*6.123233995736766e-17)).*(-t116-t202+t271+t409+t463+t474+t476+t491+t501+xg10).*-2.0-(t410+t461).*(-t116-t202+t271+t409+t463+t474+xg7).*2.0;
mt1 = [t505.*(-t116-t202+t271+t409+t463+t474+t476+t491+t501).*-2.0+t116.*t120.*2.0+(t115+t194+t269+t406+t459+t471+t472+t489+t499).*(-t116-t202+t271+t409+t463+t474+t476+t491+t501+xg10).*2.0-t115.*(t116-xg1).*2.0+t359.*(t116+t202+t273).*2.0+t502.*(t116+t202+t273-t409+t462-t474).*2.0-(t115+t194+t269).*(t116+t202+t273-xg4).*2.0+(t115+t194+t269+t406+t459+t472).*(-t116-t202+t271+t409+t463+t474+xg7).*2.0,et3+et4+et5,et6+et7+et8+et9];
mt2 = [t502.*(t196-t268+t395-t425+Ll5.*(t324-t325)).*-2.0+t276.*(t111+t211).*2.0-t505.*(t196-t268+t395+t398-t425-t494+Ll5.*(t324-t325)+Ll7.*t23.*(t393-t424)).*2.0-(t204+t270).*(t116+t202+t273-xg4).*2.0+(t204+t270+t411+t414+t436).*(-t116-t202+t271+t409+t463+t474+xg7).*2.0+t492.*(t111+t211+t350+t351+t370).*2.0-t359.*(t196-t268).*2.0+(t204+t270+t408+t411+t414+t436-t496+Ll7.*t23.*(t412-t431)).*(-t116-t202+t271+t409+t463+t474+t476+t491+t501+xg10).*2.0+t504.*(t111+t211+t350+t351+t357+t370+t444+t479).*2.0];
mt3 = [(t403+t416+t432).*(-t116-t202+t271+t409+t463+t474+xg7).*-2.0-(t403-t408+t416+t432+t483+t496).*(-t116-t202+t271+t409+t463+t474+t476+t491+t501+xg10).*2.0+t492.*(t350+t351+t370).*2.0-t502.*(t395-t425+Ll5.*(t324-t325)).*2.0-t505.*(t395+t398-t425-t494+Ll5.*(t324-t325)+Ll7.*t23.*(t393-t424)).*2.0+t504.*(t350+t351+t357+t370+t444+t479).*2.0,et10+et11+et12];
mt4 = [t505.*(Ll7.*t26.*(t455-t3.*(t324-t325))-Ll7.*t23.*t497).*2.0+(Ll7.*t23.*t500-Ll7.*t26.*(t407-t460)).*(-t116-t202+t271+t409+t463+t474+t476+t491+t501+xg10).*2.0-t504.*(Ll7.*t23.*(-t381+t382+t5.*(t250+t255).*6.123233995736766e-17)+Ll7.*t26.*t464).*2.0];
dJ_ll_dq = [mt1,mt2,mt3,mt4];
end

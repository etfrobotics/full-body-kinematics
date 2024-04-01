function H = dh_calc(theta, d, a, alpha, offset)

ct = cos(theta + offset);
st = sin(theta + offset);
ca = cos(alpha);
sa = sin(alpha);

H = [ ct    -st*ca   st*sa     a*ct ; ...
    st    ct*ca    -ct*sa    a*st ; ...
    0     sa       ca        d    ; ...
    0     0        0         1         ];
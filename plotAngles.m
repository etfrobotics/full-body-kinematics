function plotAngles(q)

l_lim = [-pi,-pi,-pi,-pi/6,-pi/2,-3*pi/4,0,-pi,-pi/4,-5*pi/6,-pi/4,-pi/4,-3*pi/4,-pi/4,-3*pi/4,-pi/6,0,-pi/6];
u_lim = [pi,pi,pi,pi/2,pi,pi/4,5*pi/6,pi/2,3*pi/4,0,3*pi/4,pi/6,0,pi/6,pi/4,pi/3,3*pi/4,pi/4];

figure('WindowState','maximized');
for ii = 1:18
    subplot(3, 6, ii)
    plot(q(:, ii), 'b');
    ylim([l_lim(ii),u_lim(ii)])
    grid on
    title(sprintf("$$q_{%d}(t)$$: %s", ii, jointIndex2Name(ii)), 'Interpreter','latex');
end
end
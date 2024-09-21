function plotAngles(t,q)

l_lim = [-pi,-pi,-pi,-pi/6,-pi/2,-3*pi/4,0,-pi,-pi/4,-5*pi/6,-pi/4,-pi/4,-3*pi/4,-pi/4,-3*pi/4,-pi/6,0,-pi/6];
u_lim = [pi,pi,pi,pi/4,pi,pi/4,5*pi/6,pi/2,3*pi/4,0,3*pi/4,pi/6,0,pi/6,pi/4,pi/3,3*pi/4,pi/4];

q(:,14) = q(:,14)+pi/6;
q(:,18) = q(:,18)-pi/6;
fig = figure('WindowState','maximized');
for ii = 1:18
    subplot(3, 6, ii)
    plot(t,q(:, ii), 'b');
    ylim([l_lim(ii),u_lim(ii)])
    grid on
    title(sprintf("$$q_{%d}(t)$$: %s", ii, jointIndex2Name(ii)), 'Interpreter','latex');
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(sprintf('Time [s]'), 'Interpreter','latex')
ylabel(sprintf('Angle [rad]'), 'Interpreter','latex')
end
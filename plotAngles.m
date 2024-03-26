function plotAngles(q)

figure('WindowState','maximized');
for ii = 1:34
    subplot(6, 6, ii)
    plot(q(:, ii), 'b');
    title(sprintf("$$q_{%d}(t)$$: %s", ii, jointIndex2Name(ii)), 'Interpreter','latex');
end
end
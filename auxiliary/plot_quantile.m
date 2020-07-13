function plot_quantile(xdata, ydata, q, color)
% plot_quantile(time, mcvals, q, color) 
% plots the quantiles of the data given in mcvals. q can be a scalar
% between 0 and 1 or a cell of scalars between 0 and 1. plot_quantile will
% plot the q-quantile as well as the (1-q)-quantile of the given data.

if ~iscell(q)
    q = {q};
end

for i = 1:length(q)
    hold on;
    xdata = xdata(:)';
    plot(xdata, quantile(ydata, q{i}/2), 'k');
    plot(xdata, quantile(ydata, 1-q{i}/2), 'k');
    fill([xdata, fliplr(xdata)], [quantile(ydata, q{i}/2), fliplr(quantile(ydata, 1-q{i}/2))], color);
    alpha(1/length(q));
end



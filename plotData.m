function plotData(xdata,ydata,titleStr,xaxis,yaxis,filename)
figure;
xdata = xdata(ydata > 0);
ydata = ydata(ydata > 0);
plot(xdata,ydata)
hold on
plot(xdata,ydata,'b.')
title(titleStr)
xlabel(xaxis)
ylabel(yaxis)
axis([0,2*pi,0, 4])
if nargin > 5
    set(gcf,'paperposition',[0 0 5 4.2])
    set(gca,'linewidth',1);
    print('-depsc','-r300',[filename '100.eps'])
end
end
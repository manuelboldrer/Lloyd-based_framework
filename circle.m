function h = circle(x,y,r,color)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;

% h = plot(xunit,yunit,'Color',color);
h = patch(xunit',yunit','r','FaceColor',color,'EdgeColor','k','FaceAlpha',0.5);
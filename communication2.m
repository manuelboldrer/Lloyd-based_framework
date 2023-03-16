function [Xvis1,Yvis1]= communication2(x,y,obs,obsD,VobsD,rr,xlim,ylim,Rmax,epsi,q,dA,Rcoh,th,b,Ad,Rmaxt0)

a=1;
b=1;
alpha=0;
trasl = 0.;
theta = -atan2(VobsD(:,1,2),VobsD(:,1,1));
Xvis = cell(1,1);
Yvis = cell(1,1);

% x(q) = x(q)+b*cos(th(q));
% y(q) = y(q)+b*sin(th(q));

%         x_egg_plot   = a*cos(0:0.1:2*pi);
%         y_egg_plot   = b*exp(-alpha*(a*cos(0:0.1:2*pi))/2).*sin(0:0.1:2*pi);
%         for qq = 1:length(x_egg_plot)
%             tmp_egg        =  [cos(theta) -sin(theta);sin(theta) cos(theta)]*[x_egg_plot(qq) + trasl; y_egg_plot(qq)];
%             x_egg_plot(qq) = tmp_egg(1);
%             y_egg_plot(qq) = tmp_egg(2);
%         end
% x_circ = zeros()
% Xvis = cell(1,1);
% Yvis = cell(1,1);
th = 0:0.07:2*pi;
p=0;
xcirc = zeros(length(th),1);
ycirc = zeros(length(th),1);

% xcirc1 = zeros(length(th),1);
% ycirc1 = zeros(length(th),1);
% for r = 0:0.05:1.5*dA
%     for j = 1:length(th)
%         
%         xcirc1(j) = x(q)+r*cos(th(j));
%         ycirc1(j) = y(q)+r*sin(th(j));
%     end
% end
for r=0:0.05:Rmaxt0
    for j = 1:length(th)
        xcirc(j) = x(q)+r*cos(th(j));
        ycirc(j) = y(q)+r*sin(th(j));
        
        if xcirc(j) >xlim(2)
            xcirc(j)=xlim(2);
        elseif xcirc(j) <= xlim(1)
            xcirc(j) = xlim(1);
        end
        if ycirc(j) >ylim(2)
            ycirc(j)=ylim(2);
        elseif ycirc(j) <= ylim(1)
            ycirc(j) = ylim(1);
        end
  if noCollision([x(q),y(q)]',[xcirc(j),ycirc(j)]',obs)==1 
            %       if noCollision([x,y]',[xcirc(j),ycirc(j)]',obs)==1 && (rr - (X-0.5)^2-(Y-0.5)^2)<0
            p=p+1;
            
            Xvis{1}(p) = xcirc(j);
            Yvis{1}(p) = ycirc(j);
            
        end
    end
end

Xvis1  = Xvis{1};
Yvis1  = Yvis{1};
end


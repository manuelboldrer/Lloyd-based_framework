function [Xvis1,Yvis1]= visibilitypoints2(x,y,ang,obs,obsD,VobsD,rr,xlim,ylim,Rmax,epsi,q,dA,Rcoh,th,b,Ad, Rmaxt0)
a=1;
b=1;
alpha=0;
trasl = 0.;
theta = -atan2(VobsD(:,1,2),VobsD(:,1,1));
Xvis = cell(1,1);
Yvis = cell(1,1);
count1 = 0;
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
th = 0:0.1:2*pi;
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
% syms x1 y1
% for jj = 1:length(x)
%         solu = solve([(x1-x(q))^2+(y1-y(q))^2-Rmax(q)^2,(x1-x(jj))^2+(y1-y(jj))^2-Rmax(jj)^2],[x1,y1]);
%                 xsol1(jj) = eval(solu.x1(1));
%                 xsol2(jj) = eval(solu.x1(2));
%                 ysol1(jj) = eval(solu.y1(1));
%                 ysol2(jj) = eval(solu.y1(2));
% %                 Xsol(jj) = .5*(xsol1+xsol2);
% %                 Ysol(jj) = .5*(ysol1+ysol2);
%
% m(jj) = (ysol2(jj)-ysol1(jj))/(xsol2(jj)-xsol1(jj));
%
% end



% d = zeros(length(x),1);
% ww = zeros(length(x),1);
% for jj = 1:length(x)
%     ww(jj) = norm([x(jj)-x(q),y(jj)-y(q)]);
%     
%     if ww(jj) > Rmax(jj)+Rmax(q)
%         d(jj) = Rmax(q)/(Rmax(jj)+Rmax(q))*ww(jj);
%     elseif  ww(jj) <= Rmax(jj)+Rmax(q) && ww(jj) >= abs(Rmax(q)-Rmax(jj))
%         d(jj) = (Rmax(q)^2-Rmax(jj)^2 + ww(jj)^2)/(2*ww(jj));
%     else
%         disp('ouch')
%     end
% end
for r=0:0.05:Rmaxt0(q)
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
        %         for w = 1:size(obsD,2)
        %         for w = 1:size(obsD,2)
        %             X(w)=xcirc(j)-obsD(:,w,1);
        %             Y(w)=ycirc(j)-obsD(:,w,2);
        %         end
        vec = zeros(size(obsD,2),1);
        %         for w = 1:size(obsD,2)
        %             vec(w,1) = obsD(:,w,3)/2-((X(w)-obsD(:,w,3)/2)*cos(theta)-(Y(w)-obsD(:,w,3)/2)*sin(theta))^2/a^2-((X(w)-obsD(:,w,3)/2)*sin(theta)+(Y(w)-obsD(:,w,3)/2)*cos(theta))^2/b^2*exp(alpha*((X(1)-obsD(:,w,3)/2)*cos(theta)-(Y(w)-obsD(:,w,3)/2)*sin(theta)))<0;
        %         end
        count = 0;
     
        vec1 = zeros(length(x)-1,1);
        vec2 = zeros(length(x)-1,1);
        
        

       
        
        for jj = 1:length(x)
            if jj ~= q
                count = count+1;
                uvec = [x(q)-x(jj);y(q)-y(jj)]/norm([x(q)-x(jj);y(q)-y(jj)]);
                m = (y(q)-y(jj))/(x(q)-x(jj));
                mr = -1/m;
                xm = .5*(x(q)+x(jj));
                ym = .5*(y(q)+y(jj));
                
                dm = norm([xm-x(q);ym-y(q)]);
%                 if dA(jj) >0.175
%                     dA(q) = Rmax(q)-0.2;
%                 else
%                     dA(q) = 0.175;
%                 end

                if dm < dA(q)+dA(jj)
                    %                     syms X Y
                    %                     vars = [X Y];
                    %                     eqs =[atan2(Y,X) == m, sqrt(X^2+Y^2) == 2*dA - dm];
                    %                     sol = solve(eqs,vars);
                    %                     if q ==3
                    %                         keyboard
                    %                     end
                    
                    sol.X = xm + (dA(q)+dA(jj)-dm)*uvec(1);
                    sol.Y = ym + (dA(q)+dA(jj)-dm)*uvec(2);
                    if norm([x(q)-x(jj),y(q)-y(jj)])>=dA(q)+dA(jj)
                        if y(q) > -1/m*(x(q)-sol.X)+sol.Y
                            vec1(count,1) = ycirc(j) > -1/m*(xcirc(j)-sol.X)+sol.Y;
                        else
                            vec1(count,1) = ycirc(j) < -1/m*(xcirc(j)-sol.X)+sol.Y;
                        end
                    else
                        if y(q) > -1/m*(x(q)-sol.X)+sol.Y
                            dm = .5*(dA(q)+dA(jj));
                            sol.X = xm + (dA(q)+dA(jj)-dm)*uvec(1);
                            sol.Y = ym + (dA(q)+dA(jj)-dm)*uvec(2);
                            vec1(count,1) = ycirc(j) < -1/m*(xcirc(j)-sol.X)+sol.Y;
                        else
                            vec1(count,1) = ycirc(j) > -1/m*(xcirc(j)-sol.X)+sol.Y;
                        end
                    end
                    %                     vec1(count,1) = ycirc(j) > -1/m*(xcirc(j)-sol.X)+sol.Y;
                    %                     plot(xcirc(j),ycirc(j),'bd')
                    
                    %                     2*dA - (xcirc(j)-x(jj))^2 - (ycirc(j)-y(jj))^2 < 0;
                else
                    vec1(count,1) = 1;
                end
%                 count1 = count1+1;
                
                %                 if y(q) > -1/m*(x(q) - Xsol(jj))+Ysol(jj)
                %                     vec2(count1,1) = ycirc(j) > -1/m*(xcirc(j)-Xsol(jj))+Ysol(jj);
                %                 else
                %                     vec2(count1,1) = ycirc(j) < -1/m*(xcirc(j)-Xsol(jj))+Ysol(jj);
                %                 end
            end
        end
        count = 0;
        for w = 1:size(obsD,2)
            count = count+1;
            uvec = [x(q)-obsD(:,w,1);y(q)-obsD(:,w,2)]/norm([x(q)-obsD(:,w,1);y(q)-obsD(:,w,2)]);
            m = (y(q)-obsD(:,w,2))/(x(q)-obsD(:,w,1));
            mr = -1/m;
            xm = .5*(x(q)+obsD(:,w,1));
            ym = .5*(y(q)+obsD(:,w,2));
            dm = norm([xm-x(q);ym-y(q)]);
            
            if dm < dA(q)+obsD(:,w,3) 
                
                sol.X = xm + (dA(q)+obsD(:,w,3)-dm)*uvec(1);
                sol.Y = ym + (dA(q)+obsD(:,w,3)-dm)*uvec(2);
                
                if y(q) > -1/m*(x(q)-sol.X)+sol.Y
                    vec(count,1) = ycirc(j) > -1/m*(xcirc(j)-sol.X)+sol.Y;
                else
                    vec(count,1) = ycirc(j) < -1/m*(xcirc(j)-sol.X)+sol.Y;
                end
                %                     vec1(count,1) = ycirc(j) > -1/m*(xcirc(j)-sol.X)+sol.Y;
                %                     plot(xcirc(j),ycirc(j),'bd')
                
                %                     2*dA - (xcirc(j)-x(jj))^2 - (ycirc(j)-y(jj))^2 < 0;
            elseif  norm([xcirc(j)-x(q),ycirc(j)-y(q)])>norm([xcirc(j)-obsD(:,w,1),ycirc(j)-obsD(:,w,2)])
                
            else
                
                vec(count,1) = 1;
            end
        end
        count2 = 0;
        for jj = 1:length(x)
            if q ~= jj         
                alpha(jj) = atan2(y(jj)-y(q),x(jj)-x(q));                
                if (xcirc(j)-x(q))*cos(alpha(jj)) + (ycirc(j)-y(q))*sin(alpha(jj)) < norm([x(q)-x(jj),y(q)-y(jj)])*Rmax(q)/(Rmax(jj)+...
                     Rmax(q)) || ( (norm([x(q)-x(jj),y(q)-y(jj)])>Rmaxt0(q)+Rmaxt0(jj) ) && norm([xcirc(j)-x(q),ycirc(j)-y(q)]) < Rmaxt0(q)) 
%                 if ((xcirc(j)-x(q))^2 + (ycirc(j)-y(q))^2) < ((xcirc(j)-x(jj))^2 + (ycirc(j)-y(jj))^2) ...
%                         || ( (norm([x(q)-x(jj),y(q)-y(jj)])>Rmaxt0(q)+Rmaxt0(jj) ) && norm([xcirc(j)-x(q),ycirc(j)-y(q)]) < Rmaxt0(q)) %Rmax(jj)/(Rmax(q)+Rmax(jj))

%                                       Rmax(q)) || ( (norm([x(q)-x(jj),y(q)-y(jj)])>Rmaxt0(q)+Rmaxt0(jj) || noCollision([x(q),y(q)]',[x(jj),y(jj)]',obs)==0 ) && norm([xcirc(j)-x(q),ycirc(j)-y(q)]) < Rmaxt0(q) ) 
% 
                    %Rmax(jj)/(Rmax(q)+Rmax(jj))
%                                if 1/Rmax(q)^2*((xcirc(j)-x(q))^2 + (ycirc(j)-y(q))^2) < 1/Rmax(jj)^2*((xcirc(j)-x(jj))^2 + (ycirc(j)-y(jj))^2) %Rmax(jj)/(Rmax(q)+Rmax(jj))
%                                if ((xcirc(j)-x(q))^2 + (ycirc(j)-y(q))^2)-Rmax(q)^2 < ((xcirc(j)-x(jj))^2 + (ycirc(j)-y(jj))^2)-Rmax(jj)^2 %Rmax(jj)/(Rmax(q)+Rmax(jj))

                    %                 if  norm([xcirc(j)-x(q),ycirc(j)-y(q)])/Rmax(q) <= norm([xcirc(j)-x(jj),ycirc(j)-y(jj)])/Rmax(jj)
                %                   if xcirc(j)-x(q)
                
                count2 = count2+1 ;
                vec2(count2,1) = 1;
                
                %             plot(xcirc(j),ycirc(j),'ro')
                %             hold on
%                 else
%                     count2 = count2+1 ;
%                     vec2(count2,1) = 1;
                end
                
%                  if  Ad(q,jj) ==0 && norm([xcirc(j)-x(q),ycirc-y(q)]) < Rmaxt0(q)
% %                     
%                      count2 = count2+1;
%                      vec2(count2,1)=1;
%                  end
            end
        end
        %% triangle
%         in = zeros(length(x),1);
%         for jj = 1:length(x)
%         
%             X1 = x(jj)+Rmax*cos(ang(jj));
%             Y1 = y(jj)+Rmax*sin(ang(jj));
%             X2 = x(jj)+Rmax*cos(ang(jj) + 120*pi/180);
%             X3 = x(jj)+Rmax*cos(ang(jj) - 120*pi/180);
%             Y2 = y(jj)+Rmax*sin(ang(jj) + 120*pi/180);
%             Y3 = y(jj)+Rmax*sin(ang(jj) - 120*pi/180);
%             X = [X1,X2,X3]';
%             Y = [Y1,Y2,Y3]';
%             kh = convhull(X,Y);
%             in(jj) = inpolygon(xcirc(j) ,ycirc(j),X,Y);
%         end
        
        %%
        coh = Rcoh - (xcirc(j)-mean(x))^2 - (ycirc(j)-mean(y))^2 > 0 ;
        
        %         if noCollision([x(q),y(q)]',[xcirc(j),ycirc(j)]',obs)==1 && vec(1) && vec(2)&& vec(3) && coh...
        %                 && vec1(1) && vec1(2) && vec1(3) && vec1(4) %&& vec1(5) && vec1(6) && vec1(7) && vec1(8) && vec1(9)%...
        %&& vec1(10) && vec1(11) && vec1(12) && vec1(13) && vec1(14) && vec1(15) && vec1(16) && vec1(17) && vec1(18) && vec1(19)
        if noCollision([x(q),y(q)]',[xcirc(j),ycirc(j)]',obs)==1 && all(vec1) && coh && all(vec2) && all(vec) 
            %       if noCollision([x,y]',[xcirc(j),ycirc(j)]',obs)==1 && (rr - (X-0.5)^2-(Y-0.5)^2)<0
            p=p+1;
            
            Xvis{1}(p) = xcirc(j);
            Yvis{1}(p) = ycirc(j);
                

            
        end
end

Xvis1  = Xvis{1};
Yvis1  = Yvis{1};
end


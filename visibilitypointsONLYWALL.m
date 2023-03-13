function [Xvis1,Yvis1]= visibilitypointsONLYWALL(x,y,obs,obsD,VobsD,rr,xlim,ylim,Rmax,epsi,q,dA,Rcoh,th,b,Ad)

a=1;
b=1;
alpha=0;
trasl = 0.;
theta = -atan2(VobsD(:,1,2),VobsD(:,1,1));
Xvis = cell(1,1);
Yvis = cell(1,1);

th = 0:0.1:2*pi;
p=0;
xcirc = zeros(length(th),1);
ycirc = zeros(length(th),1);



for r=0:0.05:Rmax(q)
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

        vec = zeros(size(obsD,2),1);

        count = 0;
        count1 = 0;
        vec1 = zeros(length(x)-1,1);
        vec2 = zeros(length(x)-1,1);
%         for jj = 1:length(x)
%             if jj ~= q
%                 count = count+1;
%                 uvec = [x(q)-x(jj);y(q)-y(jj)]/norm([x(q)-x(jj);y(q)-y(jj)]);
%                 m = (y(q)-y(jj))/(x(q)-x(jj));
%                 mr = -1/m;
%                 xm = .5*(x(q)+x(jj));
%                 ym = .5*(y(q)+y(jj));
%                 
%                 dm = norm([xm-x(q);ym-y(q)]);
%                 
%                 if dm < dA(q)+dA(jj)
%                     %                     syms X Y
%                     %                     vars = [X Y];
%                     %                     eqs =[atan2(Y,X) == m, sqrt(X^2+Y^2) == 2*dA - dm];
%                     %                     sol = solve(eqs,vars);
%                     %                     if q ==3
%                     %                         keyboard
%                     %                     end
%                     
%                     sol.X = xm + (dA(q)+dA(jj)-dm)*uvec(1);
%                     sol.Y = ym + (dA(q)+dA(jj)-dm)*uvec(2);
%                     
%                     if y(q) > -1/m*(x(q)-sol.X)+sol.Y
%                         vec1(count,1) = ycirc(j) > -1/m*(xcirc(j)-sol.X)+sol.Y;
%                     else
%                         vec1(count,1) = ycirc(j) < -1/m*(xcirc(j)-sol.X)+sol.Y;
%                     end
%                     %                     vec1(count,1) = ycirc(j) > -1/m*(xcirc(j)-sol.X)+sol.Y;
%                     %                     plot(xcirc(j),ycirc(j),'bd')
%                     
%                     %                     2*dA - (xcirc(j)-x(jj))^2 - (ycirc(j)-y(jj))^2 < 0;
%                 else
%                     vec1(count,1) = 1;
%                 end
% %                 count1 = count1+1;
%                 
%                 %                 if y(q) > -1/m*(x(q) - Xsol(jj))+Ysol(jj)
%                 %                     vec2(count1,1) = ycirc(j) > -1/m*(xcirc(j)-Xsol(jj))+Ysol(jj);
%                 %                 else
%                 %                     vec2(count1,1) = ycirc(j) < -1/m*(xcirc(j)-Xsol(jj))+Ysol(jj);
%                 %                 end
%             end
%         end
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
            else
                vec(count,1) = 1;
            end
        end
        count2 = 0;
        for jj = 1:length(x)
%             if q ~= jj          
%                 alpha(jj) = atan2(y(jj)-y(q),x(jj)-x(q));
%                 if (xcirc(j)-x(q))*cos(alpha(jj)) + (ycirc(j)-y(q))*sin(alpha(jj)) < norm([x(q)-x(jj),y(q)-y(jj)])*Rmax(q)/(Rmax(jj)+...
%                        Rmax(q)) %|| ( Ad(q,jj) ==0 && norm([xcirc(j)-x(q),ycirc(j)-y(q)]) < Rmax(q) ) 
%                
%                     %Rmax(jj)/(Rmax(q)+Rmax(jj))
% %                                if 1/Rmax(q)^2*((xcirc(j)-x(q))^2 + (ycirc(j)-y(q))^2) < 1/Rmax(jj)^2*((xcirc(j)-x(jj))^2 + (ycirc(j)-y(jj))^2) %Rmax(jj)/(Rmax(q)+Rmax(jj))
% %                                if ((xcirc(j)-x(q))^2 + (ycirc(j)-y(q))^2)-Rmax(q)^2 < ((xcirc(j)-x(jj))^2 + (ycirc(j)-y(jj))^2)-Rmax(jj)^2 %Rmax(jj)/(Rmax(q)+Rmax(jj))
% 
%                     %                 if  norm([xcirc(j)-x(q),ycirc(j)-y(q)])/Rmax(q) <= norm([xcirc(j)-x(jj),ycirc(j)-y(jj)])/Rmax(jj)
%                 %                   if xcirc(j)-x(q)
%                 count2 = count2+1 ;
%                 vec2(count2,1) = 1;
%                 %             plot(xcirc(j),ycirc(j),'ro')
%                 %             hold on
%                 end
% %             elseif Ad(q,jj) ==0
% %                 if norm([xcirc(j)-x(q),ycirc-y(q)]) < Rmax(q)
% %                     
% %                     count2 = count2+1;
% %                     vec2(count2,1)=1;
% %                 end
%             end
        
        
        coh = Rcoh - (xcirc(j)-mean(x))^2 - (ycirc(j)-mean(y))^2 > 0 ;
        
        %         if noCollision([x(q),y(q)]',[xcirc(j),ycirc(j)]',obs)==1 && vec(1) && vec(2)&& vec(3) && coh...
        %                 && vec1(1) && vec1(2) && vec1(3) && vec1(4) %&& vec1(5) && vec1(6) && vec1(7) && vec1(8) && vec1(9)%...
        %&& vec1(10) && vec1(11) && vec1(12) && vec1(13) && vec1(14) && vec1(15) && vec1(16) && vec1(17) && vec1(18) && vec1(19)
        if noCollision([x(q),y(q)]',[xcirc(j),ycirc(j)]',obs)==1  && all(vec)
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


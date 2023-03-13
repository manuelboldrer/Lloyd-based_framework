function nc = noCollision(n2, n1, o)
A = [n1(1) n1(2)];
B = [n2(1) n2(2)];
obs = zeros(size(o,1),4);
Nc  = nan(size(o,1),1);
for j = 1:size(o,1)
    obs(j,:) = [o(j,1) o(j,2) o(j,1)+o(j,3) o(j,2)+o(j,4)];
    
    C1 = [obs(j,1),obs(j,2)];
    D1 = [obs(j,1),obs(j,4)];
    C2 = [obs(j,1),obs(j,2)];
    D2 = [obs(j,3),obs(j,2)];
    C3 = [obs(j,3),obs(j,4)];
    D3 = [obs(j,3),obs(j,2)];
    C4 = [obs(j,3),obs(j,4)];
    D4 = [obs(j,1),obs(j,4)];
    
    % Check if path from n1 to n2 intersects any of the four edges of the
    % obstacle
    
    ints1 = ccw(A,C1,D1) ~= ccw(B,C1,D1) && ccw(A,B,C1) ~= ccw(A,B,D1);
    ints2 = ccw(A,C2,D2) ~= ccw(B,C2,D2) && ccw(A,B,C2) ~= ccw(A,B,D2);
    ints3 = ccw(A,C3,D3) ~= ccw(B,C3,D3) && ccw(A,B,C3) ~= ccw(A,B,D3);
    ints4 = ccw(A,C4,D4) ~= ccw(B,C4,D4) && ccw(A,B,C4) ~= ccw(A,B,D4);
    
    if ints1==0 && ints2==0 && ints3==0 && ints4==0
        Nc(j) = 0;
    else
        Nc(j) = 1;
    end
end
if sum(Nc) == 0 
    nc = 1;
else 
    nc = 0;
end
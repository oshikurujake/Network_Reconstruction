function [data1,ratio,bound_nodes,CM] = kirchhoff_grid(hn,vn,HR,VR)
% This algorithm dealing with the reconstruction process of arbitrary 
% resistor network.

% Define the horizontal and vertical elements number

% h = input('Horizontal elements = ');
% v = input('Vertical elemnts = ');

% Number of nodes in each direction
h = hn - 1;
v = vn - 1;
total_nodes = hn*vn;

% Define and sort the boundary nodes for the grid
bound_nodes = unique([1:hn,1:vn:(total_nodes-hn+1),hn:vn:total_nodes,(total_nodes-hn+1):total_nodes]);

ratio = length(bound_nodes)/total_nodes;
% Resistance to conductance
He = 1./HR;
Ve = 1./VR;

% Establish the conductance matrix
CM = zeros(hn*vn,hn*vn);
% for the point on the left top
CM(1,1) = He(1,1)+Ve(1,1);
CM(1,2) = -He(1,1);
CM(1,vn+1) = -Ve(1,1);
% for piont on the right top
CM(hn,hn) = He(1,h)+Ve(1,hn);
CM(hn,hn-1) = -He(1,h);
CM(hn,hn*2) = -Ve(1,hn);
% for piont on the left bot
CM(hn*v+1,hn*v+1) = He(vn,1)+Ve(v,1);
CM(hn*v+1, hn*(v-1)+1) = -Ve(v,1);
CM(hn*v+1, hn*v+2) = -He(vn,1);
% for piont on he right bot
CM(hn*vn,hn*vn) = He(vn,h)+Ve(v,hn);
CM(hn*vn,hn*vn-1) = -He(vn,h);
CM(hn*vn,(vn-1)*hn) = -Ve(v,hn);
% for evement hith three neighbors
% top
for n = 2:hn-1
    CM(n,n) = He(1,n-1)+He(1,n)+Ve(1,n);
    CM(n,n-1) = -He(1,n-1);
    CM(n,n+1) = -He(1,n);
    CM(n,hn+n) = - Ve(1,n);
end
% bot
for n =2:hn-1
    CM(v*hn+n,v*hn+n) = He(vn,n-1)+He(vn,n)+Ve(v,n);
    CM(v*hn+n,(v-1)*hn+n) = - Ve(v,n);
    CM(v*hn+n,v*hn+n-1) = -He(vn,n-1);
    CM(v*hn+n,v*hn+n+1) = -He(vn,n);
end
% left
for n = 2:vn-1
    CM((n-1)*hn+1,(n-1)*hn+1) = Ve(n-1,1)+Ve(n,1)+He(n,1);
    CM((n-1)*hn+1,(n-2)*hn+1) = -Ve(n-1,1);
    CM((n-1)*hn+1,n*hn+1) = -Ve(n,1);
    CM((n-1)*hn+1,(n-1)*hn+2) = -He(n,1);
end
% right
for n = 2:vn-1
    CM(n*hn,n*hn) = Ve(n-1,hn)+Ve(n,hn)+He(n,h);
    CM(n*hn,(n-1)*hn) = -Ve(n-1,hn);
    CM(n*hn,(n+1)*hn) = -Ve(n,hn);
    CM(n*hn,n*hn-1) = -He(n,h);
end

% Nodes connect hith 4 neighbors
for r = 2:v
    for c = 2:h
        CM((r-1)*hn+c,(r-1)*hn+c) = Ve(r-1,c)+Ve(r,c)+He(r,c-1)+He(r,c);
        CM((r-1)*hn+c,(r-2)*hn+c) = -Ve(r-1,c);
        CM((r-1)*hn+c,r*hn+c) = -Ve(r,c);
        CM((r-1)*hn+c,(r-1)*hn+c-1) = -He(r,c-1);
        CM((r-1)*hn+c,(r-1)*hn+c+1) = -He(r,c);
    end
end

% remove the corner nodes
CM(find(sum(abs(CM),1)==0),:)=[];
CM(:,find(sum(abs(CM),1)==0))=[];

% Calculate the potenial for each nodes accordingly.
data1 = [];
for i =1:length(bound_nodes)
    t = bound_nodes(i);
    if t == bound_nodes(length(bound_nodes))
        g = bound_nodes(1);
    else
        g = bound_nodes(i+1);
    end 
        II = zeros(total_nodes,1);
        II(t) = 1;
        II(g) = -1;
        P = zeros(total_nodes,1);
        P = CM\II;
        tv = 1;
        % Calculate the voltage
        while tv <= length(bound_nodes)
        
            if t == bound_nodes(1) 
                inj_1 = bound_nodes(length(bound_nodes));
                inj_2 = bound_nodes(1);
                inj_3 = bound_nodes(2);
            elseif t == bound_nodes(length(bound_nodes))
                inj_1 = bound_nodes(length(bound_nodes)-1);
                inj_2 = bound_nodes(length(bound_nodes));
                inj_3 = bound_nodes(1);
            else
                inj_1 = bound_nodes(i-1);
                inj_2 = bound_nodes(i);
                inj_3 = bound_nodes(i+1);
            end
                 
            if bound_nodes(tv) ~= inj_1 && bound_nodes(tv)~= inj_2&& bound_nodes(tv)~= inj_3
                if tv == length(bound_nodes)
                    data = P(bound_nodes(tv))-P(bound_nodes(1));
                    
                else
                    data = P(bound_nodes(tv))-P(bound_nodes(tv+1));
                end
                data1 = [data1;data];
                tv = tv +1;
                
            else
                tv = tv +1;
            end
        end
end

end
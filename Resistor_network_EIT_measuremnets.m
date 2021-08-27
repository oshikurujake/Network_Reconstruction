% Utilizing Kirchhoff current law to solve conductance behaviour 
%clear all
% Determine the size of the compvex nethork
h = 4;% Horizontav vength/evements number
v = 4;% Verticav vength/evements number
hn = h+1; % Horizontav nodes
vn = v+1; % Verticav nodes
nod = hn*vn;

% Benchmark resistor value
% R = 1;

% Resistors on horizontal and vertical grid
% HR = R*rand([hn,v]);
% VR = R*rand([h,vn]);

% HR = R*ones([hn,v]);
% VR = R*ones([h,vn]);

% transfer to conductance
He = 1./HR;
Ve = 1./VR;
% He(:,1) = 1e-3;
% He(:,end) = 1e-3;
% He(1,:) = 0;
% He(end,:) = 0;
% 
% Ve(1,:) = 1e-3;
% Ve(end,:) = 1e-3;
% Ve(:,1) = 0;
% Ve(:,end) = 0;
CM = zeros(hn*v,hn*v);
% Identify the boundary nodes
an = 1:vn*hn;
bn = sort([1:hn,(1:vn-1)*hn+1, (2:(vn-1))*hn, ((vn-1)*hn + 2):nod]);
% bn = [1:5,10:11,16:17,22:23,28:32];
C = cell(length(bn)*(length(bn)-1)+1,2);
C{1,1} = {'Terminal'};
C{1,2} = {'Ground'};
C{1,3} = {'Electrical Potential'};
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

% % for the point on the left top
% CM(1,1) = He(1,1)+Ve(1,1);
% CM(1,2) = -He(1,1);
% CM(1,vn+1) = -Ve(1,1);
% % for piont on the right top
% CM(vn,vn) = He(1,v)+Ve(1,vn);
% CM(vn,vn-1) = -He(1,v);
% CM(vn,vn*2) = -Ve(1,vn);
% % for piont on the veft bot
% CM(vn*h+1,vn*h+1) = He(hn,1)+Ve(h,1);
% CM(vn*h+1, vn*(h-1)+1) = -Ve(h,1);
% CM(vn*h+1, vn*h+2) = -He(hn,1);
% % for piont on he right bot
% CM(hn*vn,hn*vn) = He(hn,v)+Ve(h,vn);
% CM(hn*vn,hn*vn-1) = -He(hn,v);
% CM(hn*vn,(hn-1)*vn) = -Ve(h,vn);
% % for evement hith three neighbors
% % top
% for n = 2:vn-1
%     CM(n,n) = He(1,n-1)+He(1,n)+Ve(1,n);
%     CM(n,n-1) = -He(1,n-1);
%     CM(n,n+1) = -He(1,n);
%     CM(n,vn+n) = - Ve(1,n);
% end
% % bot
% for n =2:vn-1
%     CM(h*vn+n,h*vn+n) = He(hn,n-1)+He(hn,n)+Ve(h,n);
%     CM(h*vn+n,(h-1)*vn+n) = - Ve(h,n);
%     CM(h*vn+n,h*vn+n-1) = -He(hn,n-1);
%     CM(h*vn+n,h*vn+n+1) = -He(hn,n);
% end
% % veft
% for n = 2:hn-1
%     CM((n-1)*vn+1,(n-1)*vn+1) = Ve(n-1,1)+Ve(n,1)+He(n,1);
%     CM((n-1)*vn+1,(n-2)*vn+1) = -Ve(n-1,1);
%     CM((n-1)*vn+1,n*vn+1) = -Ve(n,1);
%     CM((n-1)*vn+1,(n-1)*vn+2) = -He(n,1);
% end
% % right
% for n = 2:hn-1
%     CM(n*vn,n*vn) = Ve(n-1,vn)+Ve(n,vn)+He(n,v);
%     CM(n*vn,(n-1)*vn) = -Ve(n-1,vn);
%     CM(n*vn,(n+1)*vn) = -Ve(n,vn);
%     CM(n*vn,n*vn-1) = -He(n,v);
% end
% 
% % Nodes connect hith 4 neighbors
% for r = 2:h
%     for c = 2:v
%         CM((r-1)*vn+c,(r-1)*vn+c) = Ve(r-1,c)+Ve(r,c)+He(r,c-1)+He(r,c);
%         CM((r-1)*vn+c,(r-2)*vn+c) = -Ve(r-1,c);
%         CM((r-1)*vn+c,r*vn+c) = -Ve(r,c);
%         CM((r-1)*vn+c,(r-1)*vn+c-1) = -He(r,c-1);
%         CM((r-1)*vn+c,(r-1)*vn+c+1) = -He(r,c);
%     end
% end

% remove the corner nodes
CM(find(sum(abs(CM),1)==0),:)=[];
CM(:,find(sum(abs(CM),1)==0))=[];
data1 = [];
for ii = 1:length(bn)
    t = bn(ii);
    if ii == 16
        
        g = bn(1);
    else
        g = bn(ii+1);
    end 
        II = zeros(nod,1);
        II(t) = 1;
        II(g) = -1;
        P = zeros(nod,1);
        P = CM\II;
        tv = 1;
        while tv <= length(bn)
            if bn(ii) ==1 
                inj_1 = 16;
                inj_2 = 1;
                inj_3 = 2;
            elseif bn(ii) == 25
                inj_1 = 24;
                inj_2 = 25;
                inj_3 = 1;
            else
                inj_1 = bn(ii-1);
                inj_2 = bn(ii);
                inj_3 = bn(ii+1);
            end
                 
            if bn(tv) ~= inj_1&& bn(tv)~= inj_2&& bn(tv)~= inj_3
                if tv == length(bn)
                    data = P(bn(tv))-P(bn(1));
                else
                    data = P(bn(tv))-P(bn(tv+1));
                end
                
                data1 = [data1;data];
                
            end
            tv = tv +1;
        end
end
% HR = HR+1 ;
% VR = VR+1 ;
% 
% % transfer to conductance
% He = 1./HR;
% Ve = 1./VR;
% 
% 
% CM = zeros(hn*v,hn*v);
% % Identify the boundary nodes
% an = 1:vn*hn;
% bn = sort([1:hn,(1:vn-1)*hn+1, (2:(vn-1))*hn, ((vn-1)*hn + 2):nod]);
% C = cell(length(bn)*(length(bn)-1)+1,2);
% C{1,1} = {'Terminal'};
% C{1,2} = {'Ground'};
% C{1,3} = {'Electrical Potential'};
% 
% 
% % for the point on the left top
% CM(1,1) = He(1,1)+Ve(1,1);
% CM(1,2) = -He(1,1);
% CM(1,vn+1) = -Ve(1,1);
% % for piont on the right top
% CM(vn,vn) = He(1,v)+Ve(1,vn);
% CM(vn,vn-1) = -He(1,v);
% CM(vn,vn*2) = -Ve(1,vn);
% % for piont on the veft bot
% CM(vn*h+1,vn*h+1) = He(hn,1)+Ve(h,1);
% CM(vn*h+1, vn*(h-1)+1) = -Ve(h,1);
% CM(vn*h+1, vn*h+2) = -He(hn,1);
% % for piont on he right bot
% CM(hn*vn,hn*vn) = He(hn,v)+Ve(h,vn);
% CM(hn*vn,hn*vn-1) = -He(hn,v);
% CM(hn*vn,(hn-1)*vn) = -Ve(h,vn);
% % for evement hith three neighbors
% % top
% for n = 2:vn-1
%     CM(n,n) = He(1,n-1)+He(1,n)+Ve(1,n);
%     CM(n,n-1) = -He(1,n-1);
%     CM(n,n+1) = -He(1,n);
%     CM(n,vn+n) = - Ve(1,n);
% end
% % bot
% for n =2:vn-1
%     CM(h*vn+n,h*vn+n) = He(hn,n-1)+He(hn,n)+Ve(h,n);
%     CM(h*vn+n,(h-1)*vn+n) = - Ve(h,n);
%     CM(h*vn+n,h*vn+n-1) = -He(hn,n-1);
%     CM(h*vn+n,h*vn+n+1) = -He(hn,n);
% end
% % veft
% for n = 2:hn-1
%     CM((n-1)*vn+1,(n-1)*vn+1) = Ve(n-1,1)+Ve(n,1)+He(n,1);
%     CM((n-1)*vn+1,(n-2)*vn+1) = -Ve(n-1,1);
%     CM((n-1)*vn+1,n*vn+1) = -Ve(n,1);
%     CM((n-1)*vn+1,(n-1)*vn+2) = -He(n,1);
% end
% % right
% for n = 2:hn-1
%     CM(n*vn,n*vn) = Ve(n-1,vn)+Ve(n,vn)+He(n,v);
%     CM(n*vn,(n-1)*vn) = -Ve(n-1,vn);
%     CM(n*vn,(n+1)*vn) = -Ve(n,vn);
%     CM(n*vn,n*vn-1) = -He(n,v);
% end
% 
% % Nodes connect hith 4 neighbors
% for r = 2:h
%     for c = 2:v
%         CM((r-1)*vn+c,(r-1)*vn+c) = Ve(r-1,c)+Ve(r,c)+He(r,c-1)+He(r,c);
%         CM((r-1)*vn+c,(r-2)*vn+c) = -Ve(r-1,c);
%         CM((r-1)*vn+c,r*vn+c) = -Ve(r,c);
%         CM((r-1)*vn+c,(r-1)*vn+c-1) = -He(r,c-1);
%         CM((r-1)*vn+c,(r-1)*vn+c+1) = -He(r,c);
%     end
% end
% data2 = [];
% for t = 1:length(bn)
%     if t == 16
%         g = bn(1);
%     else
%         g = bn(t+1);
%     end 
%         II = zeros(nod,1);
%         II(t) = 1;
%         II(g) = -1;
%         P = zeros(nod,1);
%         P = CM\II;
%         tv = 1;
%         while tv <= length(bn)
%             if tv == length(bn);
%                 data = P(bn(tv))-P(bn(1));
%             else
%                 data = P(bn(tv))-P(bn(tv+1));
%             end
%             data2 = [data2;data];
%             tv = tv +1;
%         end
%   
% end




% syms r4 r6 r7 r9 u5
% 
% eqns = [(P(5)-P(2))/r4 + (P(1)-P(2))/HR(1,1)+(P(3)-P(2))/HR(1,2) == -1,
%     (P(5)-P(6))/r7 + (P(3)-P(6))/VR(1,3) + (P(9)-P(6))/VR(2,3) ==0,
% (P(5)-P(8))/r9 + (P(9)-P(8))/HR(3,2) + (P(7)-P(8))/HR(3,1) == 0,
% (P(5)-P(4))/r6 + (P(1)-P(4))/VR(1,1) + (P(7)-P(4))/VR(2,1) == 0,
% (u5-P(2))/r4+(u5-P(6))/r7+(u5-P(8))/r9+(u5-P(4))/r6 == 0];
% 
% [A,b] = equationsToMatrix(eqns);
% eqn = [(P(5)-P(2))/r4 + (P(1)-P(2))/HR(1,1)+(P(3)-P(2))/HR(1,2) == -1,]
% S = solve(eqn)
%% Loop through all the boundaries nodes
% index = 2;
% (u5-P(2))/r4+(u5-P(6))/r7+(u5-P(8))/r9+(u5-P(4))/r6 == 0
% for i = 1%bn
%     for j = 2%bn
% 
%         if i ~= j
%             CM_1 = CM;
%             II = zeros(nod,1);
%             II(i) = 1;
%             II(j) = -1;
%             %CM_1(j,:) = [];
%             %CM_1(:,j) = [];
%             P = CM_1\II;
%             C{index,1} = i;
%             C{index,2} = j;
%             C{index,3} = P;
%             index = index + 1;
% 
%         end
%         
%     end
% end


% % Top line nodes as terminal
% for r = 1
%     % Loop ground at top nodes
%     for c = 2
%         if r == c
%             continue
%         else I = zeros(nod-1,1);
%             I(r) = 1;
%             CM_1 = CM;
%     %% Effects of terminal and ground
%     %% Remove ground points
%             CM_1(c,:) = 0;
%             CM_1(:,c) = 0;
%          CM_1(find(sum(abs(CM_1),1)==0),:)=[];
%          CM_1(:,find(sum(abs(CM_1),1)==0))=[];
% 
%     %% Assign positive value for terminal mid
%         CM_1(r,r) = 0;
%         CM_1(r,r) = -sum(CM_1(r,:));
% 
%     %% Solve the equation get potential for every point
%         P = CM_1\I;
%         Z(r,c) = P(r,1)/1;
%         end
%         
%     end
% end

% for r = vn+1:vn:(hn-1)*vn+1
%     for c = vn+1:vn:(hn-1)*vn+1
%         if r == c
%             continue
%         else I = zeros(nod-1,1);
%             I(r) = 1;
%             CM_1 = CM;
%     %% Effects of terminal and ground
%     %% Remove ground points
%         CM_1(c,:) = 0;
%         CM_1(:,c) = 0;
%         CM_1(find(sum(abs(CM_1),1)==0),:)=[];
%         CM_1(:,find(sum(abs(CM_1),1)==0))=[];
% 
%     %% Assign positive value for terminal mid
%         CM_1(r,r) = 0;
%         CM_1(r,r) = -sum(CM_1(r,:));
% 
%     %% Solve the equation get potential for every point
%         P = CM_1\I;
%         Z(r,c) = P(r,1)/1;
%         end
%         
%     end
% end


% re = real(Z);
% im = imag(Z);
% G = re./(re.^2+im.^2);
% Phase = atan(im/re)/pi*180;
% figure(1)
% loglog(omega.*abs(C).*R,abs(Z(:,:))./R./v.*(h+1)); % |Z|
% hold on
% xlabel('\omegaCR');% 
% ylabel('|Z|RL/(W+1)'); 
% 
% % figure(2)
% % loglog(omega.*abs(C).*R,G./R./v.*(h+1)); % |G|
% % hold on
% % xlabel('\omegaCR');% 
% % ylabel('G*RL/(W+1)'); 
% 
% figure(3)
% semilogx(omega.*abs(C).*R,re./R./v.*(h+1)); % real Z
% hold on
% xlabel('\omegaCR');% 
% ylabel('|real(Z)|RL/(W+1)'); 
% 
% figure(4)
% loglog(omega.*abs(C).*R,im./R./v.*(h+1)); % imag Z
% hold on
% xlabel('\omegaCR');% 
% ylabel('|r(Z)|RL/(W+1)'); 
% 
% figure(5)
% semilogx(omega,Phase); %pahse angle
% hold on
% xlabel('\omegaCR');% 
% ylabel('Phase Angle'); 



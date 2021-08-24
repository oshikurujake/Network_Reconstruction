%==========KCL============


clear all;

% Generation of a random rc nework
% =================
% =================
% 
% two random matrices to simulate the networks
% ================= 

% size of the network: element number
W=2; %  W should be even number
L=2;

Xc=0.55; % proportion of C= Xc=
seg=1; % the number of points included in the conducting lines , odd, related to W
se=floor(seg/2);
% W: height, L:length-----nodes number
SW=W+1;
SL=L+1;
% numbers of points on the curves
CP=20; 


U=zeros(SW,SL);
Nu=size(U(:));
Unn=Nu(1,1);

I=1;


for m=1:SW
    for n=1:L
        H0(m,n)=rand;
    end
end

for m=1:W
    for n=1:SL
        V0(m,n)=rand;
    end
end

for m=1:SW
    for n=1:L
        
            H1(m,n)=1/H0(m,n);  % admittance
        
    end
end

for m=1:W
    for n=1:SL
        
            V1(m,n)=1/V0(m,n);  % admittance
        
    end
end

H=H1;
V=V1;
% ==================================================
% define the source by controlling the values of elements on the boundaries
% for m=1:S
%     V(m,1)=10000000000000000000000000000000000000;
% end
% for m=1:S
%     V(m,S+1)=100000000000000000000000000000000000000;
% end


% point left top
GG(1,1)=H(1,1)+V(1,1);
GG(1,2)=-H(1,1);
GG(1,SL+1)=-V(1,1);
% point right top
GG(SL,SL)=H(1,L)+V(1,SL);
GG(SL,SL-1)=-H(1,L);
GG(SL,2*SL)=-V(1,SL);
% point left bottom
GG(W*SL+1,W*SL+1)=H(SW,1)+V(W,1);
GG(W*SL+1,W*SL+2)=-H(SW,1);
GG(W*SL+1,(W-1)*SL+1)=-V(W,1);
% point right bottom
GG(Unn,Unn)=H(SW,L)+V(W,SL);
GG(Unn,Unn-1)=-H(SW,L-1);
GG(Unn,Unn-SL)=-V(W,SL);
% top===========
for p=2:L
    %GG(2,1)=H(1,1)+H(1,2)+V(1,2);
    %GG(3,1)=H(1,2)+H(1,3)+V(1,3);
    %GG(4,1)=H(1,3)+H(1,4)+V(1,4);
    GG(p,p)=H(1,p-1)+H(1,p)+V(1,p);
    GG(p,p-1)=-H(1,p-1);
    GG(p,p+1)=-H(1,p);
    GG(p,p+SL)=-V(1,p);

end
% bottom=========
for p=2:L
    GG(Unn-SL+p,Unn-SL+p)=H(SW,p-1)+H(SW,p)+V(W,p);
    GG(Unn-SL+p,Unn-SL+p-1)=-H(SW,p-1);
    GG(Unn-SL+p,Unn-SL+p+1)=-H(SW,p);
    GG(Unn-SL+p,Unn-SL-SL+p)=-V(W,p);
end
% left===========
for p=2:W
    GG((p-1)*(SL)+1,(p-1)*(SL)+1)=V(p-1,1)+H(p,1)+V(p,1);
    GG((p-1)*(SL)+1,(p-1)*(SL)+1+1)=-H(p,1);
    GG((p-1)*(SL)+1,(p-2)*(SL)+1)=-V(p-1,1);
    GG((p-1)*(SL)+1,p*(SL)+1)=-V(p,1);
end
% right===========
for p=2:W
    GG(p*SL,p*SL)=V(p-1,SL)+H(p,L)+V(p,SL);
    GG(p*SL,p*SL-1)=-H(p,L);
    GG(p*SL,(p-1)*SL)=-V(p-1,SL);
    GG(p*SL,(p+1)*SL)=-V(p,SL);
end
% other nodes===========
for p=2:W
    for q=2:L
        GG((p-1)*(SL)+q,(p-1)*(SL)+q)=H(p,q-1)+H(p,q)+V(p-1,q)+V(p,q);
        GG((p-1)*(SL)+q,(p-1)*(SL)+q-1)=-H(p,q-1);  % H left
        GG((p-1)*(SL)+q,(p-1)*(SL)+q+1)=-H(p,q);  % H right
        GG((p-1)*(SL)+q,(p-2)*(SL)+q)=-V(p-1,q);  % V up
        GG((p-1)*(SL)+q,p*(SL)+q)=-V(p,q);  % V down
    end
end
% define the zero potential reference

%======================Boundary condition_current==================
Ix=II(1:Unn-2*seg+1,1);
Ix((floor((W)/2))*SL+1-2*se,1)=I;
% GG(:,15)=0;
%======================Boundary condition_current==================

for y=1:se
% for y =2
for o=1:Unn
    if imag(GG(o,floor((W)/2-y)*SL+1))<0 && GG(o,floor((W)/2)*SL+1)==0
        GG(o,floor((W)/2)*SL+1)=GG(o,floor((W)/2-y)*SL+1);
    elseif real(GG(o,floor((W)/2-y )*SL+1))<0 && GG(o,floor((W)/2)*SL+1)==0
        GG(o,floor((W)/2)*SL+1)=GG(o,floor((W)/2-y)*SL+1);
    elseif imag(GG(o,floor((W)/2-y)*SL+1))>0 && imag(GG(o,floor((W)/2)*SL+1))<0
        GG(o,floor((W)/2)*SL+1)=0;
        elseif imag(GG(o,floor((W)/2-y)*SL+1))>0 && real(GG(o,floor((W)/2)*SL+1))<0
        GG(o,floor((W)/2)*SL+1)=0;
    elseif real(GG(o,floor((W)/2-y)*SL+1))>0 && imag(GG(o,floor((W)/2)*SL+1))<0
        GG(o,floor((W)/2)*SL+1)=0;
        elseif imag(GG(o,floor((W)/2-y)*SL+1))>0 && real(GG(o,floor((W)/2)*SL+1))<0
        GG(o,floor((W)/2)*SL+1)=0;
    end
end
o=0;
for o=1:Unn
    if imag(GG(o,floor((W)/2+y)*SL+1))<0 && GG(o,floor((W)/2)*SL+1)==0
        GG(o,floor((W)/2)*SL+1)=GG(o,floor((W)/2+y)*SL+1);
    elseif real(GG(o,floor((W)/2+y)*SL+1))<0 && GG(o,floor((W)/2)*SL+1)==0
        GG(o,floor((W)/2)*SL+1)=GG(o,floor((W)/2+y)*SL+1);
    elseif imag(GG(o,floor((W)/2+y)*SL+1))>0 && imag(GG(o,floor((W)/2)*SL+1))<0
        GG(o,floor((W)/2)*SL+1)=0;
        elseif imag(GG(o,floor((W)/2+y)*SL+1))>0 && real(GG(o,floor((W)/2)*SL+1))<0
        GG(o,floor((W)/2)*SL+1)=0;
    elseif real(GG(o,floor((W)/2+y)*SL+1))>0 && imag(GG(o,floor((W)/2)*SL+1))<0
        GG(o,floor((W)/2)*SL+1)=0;
        elseif imag(GG(o,floor((W)/2+y)*SL+1))>0 && real(GG(o,floor((W)/2)*SL+1))<0
        GG(o,floor((W)/2)*SL+1)=0;
    end
end

end



for i=1:Unn
    if imag(GG(i,floor((W)/2)*SL+1))<0
        GG(floor((W)/2)*SL+1,i)=GG(i,floor((W)/2)*SL+1);
    elseif real(GG(i,floor((W)/2)*SL+1))<0
        GG(floor((W)/2)*SL+1,i)=GG(i,floor((W)/2)*SL+1);
    end
end


    

%======================contract the matrix size==================
% detect the ground point (the middle one) for deleting
GG(ceil((W+1)/2)*SL,:)=0;
GG(:,ceil((W+1)/2)*SL)=0;

% detect and then delete the points on the both right and left lines,
% exclusive the middle one on the left, which is the injection point.


for y=1:se

% right line

GG(ceil((W+1)/2-y)*SL,:)=0;
GG(:,ceil((W+1)/2-y)*SL)=0;
GG(ceil((W+1)/2+y)*SL,:)=0;
GG(:,ceil((W+1)/2+y)*SL)=0;
% 
% left line
GG(floor((W)/2-y)*SL+1,:)=0;
GG(:,floor((W)/2-y)*SL+1)=0;
GG(floor((W)/2+y)*SL+1,:)=0;
GG(:,floor((W)/2+y)*SL+1)=0;
end

GG(find(sum(abs(GG),1)==0),:)=[];
GG(:,find(sum(abs(GG),1)==0))=[];


% assgin the positive values for the injection point (the middle one)
GG((floor((W)/2))*SL+1-2*se,(floor((W)/2))*SL+1-2*se)=0;
GG((floor((W)/2))*SL+1-2*se,(floor((W)/2))*SL+1-2*se)=-(sum(GG(:,(floor((W)/2))*SL+1-2*se)));
%==============================================================
%==============================================================


Gx=GG;


Ux=Gx\Ix;


Z(t,1)=Ux(floor((W)/2)*SL+1-2*se,1)/I;
G(t,1)=real(Z(t,1))/(real(Z(t,1))^2+imag(Z(t,1))^2);
Ang(t,1)=atan(imag(Z(t,1))/abs(real(Z(t,1))))/pi*180;
MM(t,1)=real(Z(t,1))*WA(t,1);
MMM(t,1)=-imag(Z(t,1))*WA(t,1);
CC(t,1)=1/(abs(Z(t,1))*WA(t,1));



% 
% figure(1)
% loglog(WA(2:CP,:).*abs(C).*R,abs(CC(2:CP,:))./abs(C)./L.*(W+1)); % |Z|
% hold on
% xlabel('\omegaCR');% 
% ylabel('|Z|RL/(W+1)'); 

figure(1)
loglog(WA(2:CP,:).*abs(C).*R,abs(Z(2:CP,:))./R./L.*(W+1)); % |Z|
hold on
xlabel('\omegaCR');% 
ylabel('|Z|RL/(W+1)'); 

ccc=1/(abs(Z(CP,1))./R./L.*(W+1))/(WA(CP,1).*abs(C).*R)
% figure(2)
% semilogx(WA(2:CP,:).*abs(C).*R,Ang(2:CP,:)); % phase
% hold on
% 
% xlabel('\omegaCR');% 
% ylabel('Phase (deg)'); 
% 
% 
% figure(3)
% semilogx(WA(2:CP,:).*abs(C).*R,real(Z(2:CP,:))./R./L.*(W+1)); % real Z
% hold on
% 
% xlabel('\omegaCR');% 
% ylabel('|real(Z)|RL/(W+1)'); 
% % 
% figure(4)
% semilogx(WA(2:CP,:).*abs(C).*R,abs(imag(Z(2:CP,:)))./R./L.*(W+1)); % imag Z
% xlabel('\omegaCR');% 
% ylabel('|Z"|RL/(W+1)'); 
% hold on
% % 
% % 
% figure(5)
% loglog(WA(2:CP,:).*abs(C).*R,abs(MM(2:CP,:)./R./L.*(W+1))); % M" 
% xlabel('\omegaCR');% 
% ylabel('|M"|RL/(W+1)'); 
% hold on
% % 
% % 
% figure(6)
% plot(real(Z(2:CP,:))./R./L.*(W+1),abs(imag(Z(2:CP,:))./R./L.*(W+1))); % M" 
% xlabel('Zr');% 
% ylabel('Zi'); 
% hold on
% 
% figure(7)
% loglog(WA(2:CP,:).*abs(C).*R,abs(MMM(2:CP,:)./R./L.*(W+1))); % M" 
% xlabel('\omegaCR');% 
% ylabel('|real(M)|RL/(W+1)'); 
% hold on
% % 
% 

%         
% 
% 
% 
% 
%         
%         

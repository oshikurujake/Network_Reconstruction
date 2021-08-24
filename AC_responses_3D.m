
%% KCL initiation
% =================
% Import of the data.
dir ='/Users/chonpguzhai/Google_Drive/Usyd/AC/AC_code';
cd(dir)
load('Conf.mat')
A=conf;

for cyc=1 % repeat the process for * times.
% ----------------------------------------------------------------------------------------------------------------------------------------
% There are three columns in "con", clomn 1 ?particle index? and 2 (particle index) shows the all the contacts
% thorughout this random acking (or other pakcing you have), obtained from DEM. 
% In this test file, "con", there are 622 particles from this DEM simulaiton with particle index 1 to 622 
% Index -3, and index -4 in column  2 are top and bottom boundaries. 
% -1000 is the side boundaries. In this simulation, these side boundaires
% is assumed to to ideal isolator, thus being ignored.

% ----------------------------------------------------------------------------------------------------------------------------------------
A(:,4) = rand(numel(conf(:,1)),1);
A(:,5) = rand(numel(conf(:,1)),1);
% Colum 4 and 5 were randomly generated in the range of [0,1].
% Columns 4 and 5 are used for R and C values, respectivley, as they are dtermined(here the R or C elements are randomly by the contact force, column 5).
% ----------------------------------------------------------------------------------------------------------------------------------------
% column 3 is the contact force, realised from DEM. A exponential decay can
% be found for contact force distribtuions.
% Other distribtuiosn such as Guassin distribution,

% computed by,
% A(:,3) = normrnd(1,0.3,[numel(con(:,1)),1]); 

% though for granular materails, normal distribution doesn't sound reasonable.
% Also, R, C, RC, R-R//C or C-R//C should be dtermiend based on the contact properties, force and contact law, surface strucutures.
% ----------------------------------------------------------------------------------------------------------------------------------------

% ----------------------------------------------------------------------------------------------------------------------------------------
% find the number of particles in the packing.
Num_particle=max(A(:,2)); 
Fc=2.5; % this Fc is the critical force value, thresholding the contact force (colunm 5), to determine the contact to C or R type. 
% In this test, we assume the contact is either R or C.
% In the appendix of EML(2018) paper of our group, different R//c, R-R//C,
% C-R//C were considerred, with random values, but not based on the contact
% force, contact law or surface strucutres.
% ----------------------------------------------------------------------------------------------------------------------------------------
C=1E-9i; % 1E-9i capacitance (unit) magnitude, nF
R=1000; % 1000  resistance (unit) magnitude, kohm
I=1; % assume I=1A, injected on the top of the container
CP=30; % numbers of fequency points on a single curve
Conn=size(A,1); % Contact number throughout the pakcing.
Ux=zeros(FEn,CP);  % for saving the votalge

% ----------------------------------------------------------------------------------------------------------------------------------------
% Build the matrix (for all contacts including top and bottom boundaries)
% Kn records all contacts and use a the contact force (column5 of con) to represent this contact, in order to dtermine the element type
% This creteria can be changed by incorporating surface strucutres.
Kn=zeros(Num_particle+2,Num_particle+2); 
for m=1:1:Conn
    if A(m,2)>0 % if a contact is between two particles
        Kn(A(m,1),A(m,2))=A(m,3);
        Kn(A(m,2),A(m,1))=Kn(A(m,1),A(m,2));
    elseif A(m,2)==-3
        Kn(En+1,A(m,1))=A(m,3);
        Kn(A(m,1),En+1)= Kn(En+1,A(m,1));% top electrode
    elseif A(m,2)==-4
        Kn(En+2,A(m,1))=A(m,3);  % bottom electrode with the voltage being 0
        Kn(A(m,1),En+2)=Kn(En+2,A(m,1));
    end
end
% Smiliarity to Kn, Knr records all contacts and use a the R factor
% (column3 of con) to represent this contact.
Knr=zeros(Num_particle+2,Num_particle+2); 
for m=1:1:Conn
    if A(m,2)>0 % if a contact is between two particles
        Knr(A(m,1),A(m,2))=A(m,4);
        Knr(A(m,2),A(m,1))=Knr(A(m,1),A(m,2));
    elseif A(m,2)==-3
        Knr(En+1,A(m,1))=A(m,4);
        Knr(A(m,1),En+1)= Knr(En+1,A(m,1));% top electrode
    elseif A(m,2)==-4
        Knr(En+2,A(m,1))=A(m,4);  % bottom electrode with the voltage being 0
        Knr(A(m,1),En+2)=Knr(En+2,A(m,1));
    end
end
% Smiliarity to Kn, Knr records all contacts and use a the C factor
% (column3 of con) to represent this contact.
Knc=zeros(Num_particle+2,Num_particle+2); 
for m=1:1:Conn
    if A(m,2)>0 % if a contact is between two particles
        Knc(A(m,1),A(m,2))=A(m,5);
        Knc(A(m,2),A(m,1))=Knc(A(m,1),A(m,2));
    elseif A(m,2)==-3
        Knc(En+1,A(m,1))=A(m,5);
        Knc(A(m,1),En+1)= Knc(En+1,A(m,1));% top electrode
    elseif A(m,2)==-4
        Knc(En+2,A(m,1))=A(m,5);  % bottom electrode with the voltage being 0
        Knc(A(m,1),En+2)=Knc(En+2,A(m,1));
    end
end

% impedance 
Z=zeros(CP,1);
% frequency
W=zeros(CP,1);
% starting frequency
W(1,1)=1; 
% conductance
G=zeros(CP,1);
Ang=zeros(CP,1);
% frequency-dependent response;
% ----------------------------------------------------------------------------------------------------------------------------------------
for t=2:CP
    W(t,1)=W(t-1,1)*2; % frequency values in crease by a factor of 2 in order.
    GG=Kn;% GG conver the contact network, Kn, to conduction network by adding C and R values
for m=1:1:En+2
    for n=1:1:En+2
        if GG(m,n)>Fc
            GG(m,n)=-1/R*Knr(m,n);% put the R values to condcuction network, be careful of the sign
        elseif (GG(m,n)<Fc) % &(GG(m,n)>0)
            GG(m,n)=-(W(t,1)*C)*Knc(m,n);% put the C values to condcuction network, be careful of the sign
        end
    end
end
% define diagonal elements, based on the definition of KCL matrix,which is essential 
% for current balance of current input and outpu a node.
for m=1:1:En+2
    GG(m,m)=0;
    a=-sum(GG(m,:));
    GG(m,m)=a;
end
% ================== delete those nodes having zero contact=============
fr=find(sum(abs(GG),1)==0);
GG(find(sum(abs(GG),1)==0),:)=[];
GG(:,find(sum(abs(GG),1)==0))=[];
% define the zero potential reference
FEN=size(GG(:,1)); % final size of the network used
FEn=FEN(1,1)-1; % taking out groud, i.e.,-4, where all potential vlaues of nodes (particle) cotnact with bottom (-4) is defined to be 0.
Gx=zeros(FEn,FEn);   % final conduciton matrix
Ix=zeros(FEn,1);   % for saving current 
%======================Boundary condition_current==================
Ix(FEn,1)=I;% the current is injected from the top electrode boundary, i.e., -3
Gx=GG(1:FEn,1:FEn);
%======================Impedance and conduction==================
ux=Gx\Ix;
Ux(:,t)=ux;
Z(t,1)=ux(FEn,1)/I;
G(t,1)=real(Z(t,1))/(real(Z(t,1))^2+imag(Z(t,1))^2); % conductance
Ang(t,1)=atan(imag(Z(t,1))/abs(real(Z(t,1))))/pi*180; % impedance angle
% with conductance and impedance angle, you could then get other parameters
% of interest, capacitance, Q index, real or imaginary. 
end
 
figure(1)
loglog(W(2:CP,1).*abs(C).*R,abs(Z(2:CP,1))./R,'.-','linewidth',2,'markersize',10,'markerfacecolor','w'); % |Z|
hold on
xl=xlabel('$\omega^*$','interpreter','latex');% 
yl=ylabel('$|Z|/R$','interpreter','latex'); 
set([xl,yl,gca],'fontsize',24,'linewidth',2);
pbaspect([1.2,1,1]);
end




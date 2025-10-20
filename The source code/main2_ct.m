
clear;
clc;
close all;

run network_attack.m
% G = graph(link); 

edgenum = nnz(link)/2;

nodenum = size(link,1);

degree = sum(link,2);
degree1=1/max(degree);

D=diag(degree); % degree of each node
La=D-link;


T=0.3;
w1=0.4;
A=[1 sin(w1*T)/w1 0 -(1-cos(w1*T))/w1;0 cos(w1*T) 0 -sin(w1*T);0 (1-cos(w1*T))/w1 1 sin(w1*T)/w1;0 sin(w1*T) 0 cos(w1*T)]; 
G=[T^2/2 0;T 0;0 T^2/2;0 T];
L=100;    
ebu=0.02; 

x0=[60; 4.5; 40; 5.6];% m; m/s; m; m/s
xe0=[58; 5.8; 40; 3.8];

P0=2*eye(4); 

LR=[40 50 60 600];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:L+1
Q=0.1*eye(2); 
w{k}=mvnrnd([0,0],Q)';

if k==1
x{k}=x0;
else  
x{k}=A*x{k-1}+G*w{k-1};
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distance1=zeros(nodenum,L);
dd=length(LR);

MC=100; 
R1_P=15; 

MSE1_pos=zeros(L+1,MC,dd);
MSE1_vel=zeros(L+1,MC,dd);

for s=1:dd
for k=1:L
for i=1:nodenum 
R1{i}=eye(2);
v{i,k,s}=mvnrnd([0,0],R1{i})';
H{i}=[1 0 0 0;0 0 1 0];
y{i,k,s}=H{i}*x{k}+v{i,k,s};
r{i,s}=LR(s);
end

for i=1:nodenum 
distance1(i,k)=norm([x{k}(1) x{k}(3)]-[netx(i,1) nety(i,1)]);
  if distance1(i,k)<=r{i,s}
lamda{i,k,s}=1;
  else
lamda{i,k,s}=0;  
  end
end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MSE_P=zeros(L,dd);
MSE_P1=zeros(L,dd);

for s=1:dd
for Mont=1:MC
for k=1:L

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if k==1
for i=1:nodenum 
gama{i,k,s}=binornd(1,1);
end
    else
for i=1:nodenum
act=0.9;    
gama{i,k,s}=binornd(1,act);
end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dosProb=0.2;
L1{k,s}=binornd(1,dosProb,nodenum).*link;
L2{k,s}=triu(L1{k,s})+triu(L1{k,s})'; 
for i=1:nodenum
    for j=1:nodenum
        if L2{k,s}(i,j)==1
            phy{i,j,k,s}=1;
        else
            phy{i,j,k,s}=0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rou=6;                                                                                                                                  
fdiProb=0.3;
L3{k,s}=binornd(1,fdiProb,nodenum).*link;
without_L3{k,s}=link-L3{k,s};

for i=1:nodenum
    for j=1:nodenum
M{i,j}=eye(4); 
M1{i,j}=eye(2);  
aerfa{i,j,k,s}=mvnrnd([0,0,0,0],M{i,j})'; 
aerfa1{i,j,k,s}=mvnrnd([0,0],M1{i,j})';
    end
end

rand1{k,s}=rand;

if k==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Attack with detection
for i=1:nodenum
    for j=1:nodenum
  if L3{k,s}(i,j)==1 && i~=j && lamda{j,k,s}==1 && gama{j,k,s}==1
  z{j,k,s}=y{j,k,s}-H{j}*xe0;
  send_xe{i,j,k,s}=xe0+rou*aerfa{i,j,k,s};
  send_z{i,j,k,s}=z{j,k,s}+rou*aerfa1{i,j,k,s};
  Phy{i,j,k,s}=1-exp(-1/2*(send_z{i,j,k,s})'*inv(R1_P)*send_z{i,j,k,s});
     if Phy{i,j,k,s}>=rand1{k,s}
     eta{i,j,k,s}=0;
     else
     eta{i,j,k,s}=1;
     end

  elseif L3{k,s}(i,j)==1 && i~=j && lamda{j,k,s}==0 && gama{j,k,s}==1
  send_xe{i,j,k,s}=xe0+rou*aerfa{j,k}; 
  det_z{i,j,k,s}=lamda{i,k,s}*y{i,k,s}-H{i}*send_xe{i,j,k,s};
  Phy{i,j,k,s}=1-exp(-1/2 *(det_z{i,j,k,s})'*inv(R1_P)*det_z{i,j,k,s});
     if Phy{i,j,k,s}>=rand1{k,s}
     eta{i,j,k,s}=0;
     else
     eta{i,j,k,s}=1;
     end

  elseif without_L3{k,s}(i,j)==1 && i~=j && lamda{j,k,s}==1 && gama{j,k,s}==1
  q{i,j,k,s}=0;
  z{j,k,s}=y{j,k,s}-H{j}*xe0;
  send_xe{i,j,k,s}=xe0;
  send_z{i,j,k,s}=z{j,k,s};
  Phy{i,j,k,s}=1-exp(-1/2 *(send_z{i,j,k,s})'*inv(R1_P)*send_z{i,j,k,s});
     if Phy{i,j,k,s}>=rand1{k,s}
     eta{i,j,k,s}=0;
     else
     eta{i,j,k,s}=1;
     end

 elseif without_L3{k,s}(i,j)==1 && i~=j && lamda{j,k,s}==0 && gama{j,k,s}==1
  send_xe{i,j,k,s}=xe0;
  det_z{i,j,k,s}=lamda{i,k,s}*y{i,k,s}+(1-lamda{i,k,s})*H{i}*xe0-H{i}*send_xe{i,j,k,s};
  Phy{i,j,k,s}=1-exp(-1/2 *(det_z{i,j,k,s})'*inv(R1_P)*det_z{i,j,k,s});
     if Phy{i,j,k,s}>=rand1{k,s}
     eta{i,j,k,s}=0;
     else
     eta{i,j,k,s}=1;
     end
  
  else
  send_xe{i,j,k,s}=0;
  eta{i,j,k,s}=0;
  end
  end
  end

for i=1:nodenum
    att_suma=0;
    for j=1:nodenum
        if  send_xe{i,j,k,s}~=0 & eta{i,j,k,s}==1
    att_suma1=ebu*A*(1-phy{i,j,k,s})*(xe0-send_xe{i,j,k,s});
        else
    att_suma1=0;
        end
    att_suma=att_suma+att_suma1;
    end
    attack_suma{i,k,s}=att_suma;
    clear att_suma
end

delta_t=0.2;
for i=1:nodenum
P{i,k,s}=P0;
xe{i,k,s}=xe0;
yl{i,k,s}=lamda{i,k,s}*y{i,k,s}+(1-lamda{i,k,s})*H{i}*xe{i,k,s};
K{i,k,s}=(1+delta_t)*lamda{i,k,s}*A*P{i,k,s}*H{i}'*inv((1+delta_t)*H{i}*P{i,k,s}*H{i}'+R1{i});
xe{i,k+1,s}=A*xe{i,k,s}+K{i,k}*(yl{i,k,s}-H{i}*xe{i,k,s})-attack_suma{i,k,s};
F{i,k,s}=attack_suma{i,k,s};
P{i,k+1,s}=(1+delta_t)*(A-lamda{i,k,s}*K{i,k,s}*H{i})*P{i,k,s}*(A-lamda{i,k,s}*K{i,k,s}*H{i})'...
+lamda{i,k,s}*K{i,k,s}*R1{i}*K{i,k,s}'+(1+1/delta_t)*F{i,k,s}*F{i,k,s}'+G*Q*G';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Attack with detection
for i=1:nodenum
    for j=1:nodenum
  if L3{k,s}(i,j)==1 && i~=j && lamda{j,k,s}==1 && gama{j,k,s}==1
  z{j,k,s}=y{j,k,s}-H{j}*xe{j,k,s};
  send_xe{i,j,k,s}=xe{j,k,s}+rou*aerfa{i,j,k,s};
  send_z{i,j,k,s}=z{j,k,s}+rou*aerfa1{i,j,k,s};
  Phy{i,j,k,s}=1-exp(-1/2 *(send_z{i,j,k,s})'*inv(R1_P)*send_z{i,j,k,s});
     if Phy{i,j,k,s}>=rand1{k,s}
     eta{i,j,k,s}=0;
     else
     eta{i,j,k,s}=1;
     end

  elseif L3{k,s}(i,j)==1 && i~=j && lamda{j,k,s}==0 && gama{j,k,s}==1
  send_xe{i,j,k,s}=xe{j,k,s}+rou*aerfa{i,j,k,s}; 
  det_z{i,j,k,s}=lamda{i,k,s}*y{i,k,s}-H{i}*send_xe{i,j,k,s};
  Phy{i,j,k,s}=1-exp(-1/2 *(det_z{i,j,k,s})'*inv(R1_P)*det_z{i,j,k,s});
     if Phy{i,j,k,s}>=rand1{k,s}
     eta{i,j,k,s}=0;
     else
     eta{i,j,k,s}=1;
     end

 elseif without_L3{k,s}(i,j)==1 && i~=j && lamda{j,k,s}==1 && gama{j,k,s}==1
  z{j,k,s}=y{j,k,s}-H{j}*xe{j,k,s};
  send_xe{i,j,k,s}=xe{j,k,s};
  send_z{i,j,k,s}=z{j,k,s};
  Phy{i,j,k,s}=1-exp(-1/2 *(send_z{i,j,k,s})'*inv(R1_P)*send_z{i,j,k,s});
     if Phy{i,j,k,s}>=rand1{k,s}
     eta{i,j,k,s}=0;
     else
     eta{i,j,k,s}=1;
     end

 elseif without_L3{k,s}(i,j)==1 && i~=j && lamda{j,k,s}==0 && gama{j,k,s}==1
  send_xe{i,j,k,s}=xe{j,k,s};
  det_z{i,j,k,s}=lamda{i,k,s}*y{i,k,s}-H{i}*send_xe{i,j,k,s};
  Phy{i,j,k,s}=1-exp(-1/2 *(det_z{i,j,k,s})'*inv(R1_P)*det_z{i,j,k,s});
     if Phy{i,j,k,s}>=rand1{k,s}
     eta{i,j,k,s}=0;
     else
     eta{i,j,k,s}=1;
     end

  else
  send_xe{i,j,k,s}=0;
  eta{i,j,k,s}=0;
  end
    end
end

for i=1:nodenum
    att_suma=0;
    for j=1:nodenum
        if  send_xe{i,j,k,s}~=0 & eta{i,j,k,s}==1
    att_suma1=ebu*A*(1-phy{i,j,k,s})*(xe{i,k,s}-send_xe{i,j,k,s});
        else
    att_suma1=0;
        end
    att_suma=att_suma+att_suma1;
    end
    attack_suma{i,k,s}=att_suma;
    clear att_suma
end

delta_t=0.2;
for i=1:nodenum
yl{i,k,s}=lamda{i,k,s}*y{i,k,s}+(1-lamda{i,k,s})*H{i}*xe{i,k,s};
K{i,k,s}=(1+delta_t)*lamda{i,k,s}*A*P{i,k,s}*H{i}'*inv((1+delta_t)*H{i}*P{i,k,s}*H{i}'+R1{i});
xe{i,k+1,s}=A*xe{i,k,s}+K{i,k,s}*(yl{i,k,s}-H{i}*xe{i,k,s})-attack_suma{i,k,s};
F{i,k,s}=attack_suma{i,k,s};
P{i,k+1,s}=(1+delta_t)*(A-lamda{i,k,s}*K{i,k,s}*H{i})*P{i,k,s}*(A-lamda{i,k,s}*K{i,k,s}*H{i})'...
+lamda{i,k,s}*K{i,k,s}*R1{i}*K{i,k,s}'+(1+1/delta_t)*F{i,k,s}*F{i,k,s}'+G*Q*G';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MSE
x_value=cell2mat(x);

for i=1:L+1
suma1=0;

for j=1:nodenum
suma1=suma1+(xe{j,i,s}-x{i}).*(xe{j,i,s}-x{i});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MSE1_pos(i,Mont,s)=suma1(1,:)+suma1(3,:); 
MSE1_vel(i,Mont,s)=suma1(2,:)+suma1(4,:); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end

for o=1:L
sumb=0;

for MC=1:Mont
sumb=sumb+MSE1_pos(o,Mont,s);
end
    MSE_P(o,s)=sumb/MC/nodenum;
end

for o=1:L
sumb1=0;

for Mont=1:MC
sumb1=sumb1+MSE1_vel(o,Mont,s);

end
    MSE_P1(o,s)=sumb1/MC/nodenum;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
plot(1:L, MSE_P(:,1),'k-x',1:L,MSE_P(:,2),'b-.',1:L,MSE_P(:,3),'r-',1:L,MSE_P(:,4),'c-.','linewidth',1)
legend({'Sengsing range 40','Sengsing range 50','Sengsing range 60','Without limitation'},'FontName','Arial');
xlabel('Time')
ylabel('MSE^{pos}_k')
xlim([0, 100])
set(gca,'FontName','Arial');

figure(5)
plot(1:L,MSE_P1(:,1),'k-x',1:L,MSE_P1(:,2),'b-.',1:L,MSE_P1(:,3),'r-',1:L,MSE_P1(:,4),'c-.','linewidth',1)
legend({'Sengsing range 40','Sengsing range 50','Sengsing range 60','Without limitation'},'FontName','Arial');
xlabel('Time')
ylabel('MSE^{vel}_k')
xlim([0, 100])
set(gca,'FontName','Arial');
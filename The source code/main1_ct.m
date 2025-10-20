
clear;
clc;
close all;

run network_attack.m
% G = graph(link); 

edgenum = nnz(link)/2;

nodenum = size(link,1);

degree = sum(link,2);
degree1=1/max(degree);

D = diag(degree); % degree of each node
La=D-link;

T=0.3;
w1=0.4;
A=[1 sin(w1*T)/w1 0 -(1-cos(w1*T))/w1;0 cos(w1*T) 0 -sin(w1*T);0 (1-cos(w1*T))/w1 1 sin(w1*T)/w1;0 sin(w1*T) 0 cos(w1*T)]; 
G=[T^2/2 0;T 0;0 T^2/2;0 T];
L=100;    
ebu=0.02; 
delta=0.2;

x0=[60; 4.5; 40; 5.6];% m; m/s; m; m/s
xe0=[58; 5.5; 46; 4.8];

P0=2*eye(4); 


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

for k=1:L
for i=1:nodenum 
R1{i}=eye(2);
v{i,k}=mvnrnd([0,0],R1{i})';
H{i}=[1 0 0 0;0 0 1 0];
y{i,k}=H{i}*x{k}+v{i,k};
r{i}=70;
end

for i=1:nodenum 
distance1(i,k)=norm([x{k}(1) x{k}(3)]-[netx(i,1) nety(i,1)]);
if distance1(i,k)<=r{i}
lamda{i,k}=1;
else
lamda{i,k}=0;  
end
end
end

LP=40;
FP=zeros(1,LP);
FN=zeros(1,LP);
FPa=zeros(1,LP);
FNa=zeros(1,LP);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MC=500;
MSE1_pos=zeros(1,L+1,MC);
MSE2_pos=zeros(1,L+1,MC);
MSE3_pos=zeros(1,L+1,MC);
MSEa_pos=zeros(1,L+1,MC);
MSEb_pos=zeros(1,L+1,MC);

MSE1_vel=zeros(1,L+1,MC);
MSE2_vel=zeros(1,L+1,MC);
MSE3_vel=zeros(1,L+1,MC);
MSEa_vel=zeros(1,L+1,MC);
MSEb_vel=zeros(1,L+1,MC);

for Mont=1:MC
for R1_P=15
for k=1:L
    if k==1
for i=1:nodenum 
gama{i,k}=binornd(1,1);
end
    else
for i=1:nodenum
act=0.9;    
gama{i,k}=binornd(1,act);
end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Dos attack
dosProb=0.2;
L1{k}=binornd(1,dosProb,nodenum).*link;
L2{k}=triu(L1{k})+triu(L1{k})'; %发生Dos攻击的边
for i=1:nodenum
    for j=1:nodenum
        if L2{k}(i,j)==1
            phy{i,j,k}=1;
        else
            phy{i,j,k}=0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FDI attack
rou=6;                                                                                                                                  
fdiProb=0.3;
L3{k}=binornd(1,fdiProb,nodenum).*link; %发生fdi攻击的边
without_L3{k}=link-L3{k};

for i=1:nodenum
    for j=1:nodenum
M{i,j}=eye(4); 
M1{i,j}=eye(2);  
aerfa{i,j,k}=mvnrnd([0,0,0,0],M{i,j})'; 
aerfa1{i,j,k}=mvnrnd([0,0],M1{i,j})';
    end
end

rand1{k}=rand;

if k==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Attack with detection
for i=1:nodenum
    for j=1:nodenum
  if L3{k}(i,j)==1 && i~=j && lamda{j,k}==1 && gama{j,k}==1
  q{i,j,k}=1;
  z{j,k}=y{j,k}-H{j}*xe0;
  send_xe{i,j,k}=xe0+rou*q{i,j,k}*aerfa{i,j,k};
  send_z{i,j,k}=z{j,k}+rou*q{i,j,k}*aerfa1{i,j,k};
  Phy{i,j,k}=1-exp(-1/2*(send_z{i,j,k})'*inv(R1_P)*send_z{i,j,k});
     if Phy{i,j,k}>=rand1{k}
     eta{i,j,k}=0;
     else
     eta{i,j,k}=1;
     end

  elseif L3{k}(i,j)==1 && i~=j && lamda{j,k}==0 && gama{j,k}==1
  q{i,j,k}=1;
  send_xe{i,j,k}=xe0+rou*q{i,j,k}*aerfa{j,k}; 
  det_z{i,j,k}=lamda{i,k}*y{i,k}+(1-lamda{i,k})*H{i}*xe0-H{i}*send_xe{i,j,k};
  Phy{i,j,k}=1-exp(-1/2 *(det_z{i,j,k})'*inv(R1_P)*det_z{i,j,k});
     if Phy{i,j,k}>=rand1{k}
     eta{i,j,k}=0;
     else
     eta{i,j,k}=1;
     end

  elseif without_L3{k}(i,j)==1 && i~=j && lamda{j,k}==1 && gama{j,k}==1
  q{i,j,k}=0;
  z{j,k}=lamda{j,k}*(y{j,k}-H{j}*xe0);
  send_xe{i,j,k}=xe0;
  send_z{i,j,k}=z{j,k};
  Phy{i,j,k}=1-exp(-1/2 *(send_z{i,j,k})'*inv(R1_P)*send_z{i,j,k});
     if Phy{i,j,k}>=rand1{k}
     eta{i,j,k}=0;
     else
     eta{i,j,k}=1;
     end

 elseif without_L3{k}(i,j)==1 && i~=j && lamda{j,k}==0 && gama{j,k}==1
  q{i,j,k}=0;
  send_xe{i,j,k}=xe0;
  det_z{i,j,k}=lamda{i,k}*y{i,k}+(1-lamda{i,k})*H{i}*xe0-H{i}*send_xe{i,j,k};
  Phy{i,j,k}=1-exp(-1/2 *(det_z{i,j,k})'*inv(R1_P)*det_z{i,j,k});
     if Phy{i,j,k}>=rand1{k}
     eta{i,j,k}=0;
     else
     eta{i,j,k}=1;
     end
  
  else
  q{i,j,k}=0;
  send_xe{i,j,k}=0;
  eta{i,j,k}=0;
  
  end
  end
  end


for i=1:nodenum
    att_suma=0;
    for j=1:nodenum
        if  send_xe{i,j,k}~=0 & eta{i,j,k}==1
    att_suma1=ebu*A*(1-phy{i,j,k})*(xe0-send_xe{i,j,k});
        else
    att_suma1=0;
        end
    att_suma=att_suma+att_suma1;
    end
    attack_suma{i,k}=att_suma;
    clear att_suma
end

for i=1:nodenum
P{i,k}=P0;
xe{i,k}=xe0;
yl{i,k}=lamda{i,k}*y{i,k}+(1-lamda{i,k})*H{i}*xe{i,k};
K{i,k}=(1+delta)*A*P{i,k}*H{i}'*inv((1+delta)*H{i}*P{i,k}*H{i}'+R1{i});
xe{i,k+1}=A*xe{i,k}+K{i,k}*(yl{i,k}-H{i}*xe{i,k})-attack_suma{i,k};
F{i,k}=attack_suma{i,k};
P{i,k+1}=(1+delta)*(A-lamda{i,k}*K{i,k}*H{i})*P{i,k}*(A-lamda{i,k}*K{i,k}*H{i})'+lamda{i,k}*K{i,k}*R1{i}*K{i,k}'+(1+1/delta)*F{i,k}*F{i,k}'+G*Q*G';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%No attack without detection
for i=1:nodenum
  for j=1:nodenum
  if link(i,j)==1 && gama{j,k}==1
  send_xes{i,j,k}=xe0;
  else
  send_xes{i,j,k}=0;
  end
  end
end

for i=1:nodenum
    att_sum_safe=0;
    for j=1:nodenum 
        if send_xes{i,j,k}~=0
    att_sum_safe1=ebu*A*(xe0-send_xes{i,j,k});
        else
    att_sum_safe1=0;
        end
    att_sum_safe=att_sum_safe+att_sum_safe1;
    end
    attack_sums{i,k}=att_sum_safe;
    clear att_sum_safe
end
for i=1:nodenum
P_safe{i,k}=P0;
xe_safe{i,k}=xe0;
yl_safe{i,k}=lamda{i,k}*y{i,k}+(1-lamda{i,k})*H{i}*xe_safe{i,k};
K_safe{i,k}=(1+delta)*A*P_safe{i,k}*H{i}'*inv((1+delta)*H{i}*P_safe{i,k}*H{i}'+R1{i});
xe_safe{i,k+1}=A*xe_safe{i,k}+K_safe{i,k}*(yl_safe{i,k}-H{i}*xe_safe{i,k})-attack_sums{i,k};
P_safe{i,k+1}=(1+delta)*(A-lamda{i,k}*K_safe{i,k}*H{i})*P_safe{i,k}*(A-lamda{i,k}*K_safe{i,k}*H{i})'+lamda{i,k}*K_safe{i,k}*R1{i}*K_safe{i,k}'+(1+1/delta)*attack_sums{i,k}*attack_sums{i,k}'+G*Q*G';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%No attack with detection
for i=1:nodenum
    for j=1:nodenum
  if link(i,j)==1 && i~=j && lamda{j,k}==1 && gama{j,k}==1
  z1{j,k}=y{j,k}-H{j}*xe0;
  send_xe1{i,j,k}=xe0;
  send_z1{i,j,k}=z1{j,k};
  Phy1{i,j,k}=1-exp(-1/2 *send_z1{i,j,k}'*inv(R1_P)*send_z1{i,j,k});
     if Phy1{i,j,k}>=rand1{k}
     eta1{i,j,k}=0;
     else
     eta1{i,j,k}=1;
     end

  elseif link(i,j)==1 && i~=j && lamda{j,k}==0 && gama{j,k}==1
  z1{j,k}=0;
  send_xe1{i,j,k}=xe0; 
  det_z1{i,j,k}=lamda{i,k}*y{i,k}-H{i}*send_xe1{i,j,k};
  Phy1{i,j,k}=1-exp(-1/2* (det_z1{i,j,k})'*inv(R1_P)*det_z1{i,j,k});
     if Phy1{i,j,k}>=rand1{k}
     eta1{i,j,k}=0;
     else
     eta1{i,j,k}=1;
     end
  
  else
  send_xe1{i,j,k}=0;
  eta1{i,j,k}=0;

  end
  end
end

for i=1:nodenum
    att_sum_witha=0;
    for j=1:nodenum
        if send_xe1{i,j,k}~=0 & eta1{i,j,k}==1
    att_sum_witha1=ebu*A*(xe0-send_xe1{i,j,k});
        else
    att_sum_witha1=0;
        end
        att_sum_witha=att_sum_witha+att_sum_witha1;
    end
    attack_sum_witha{i,k}=att_sum_witha;
    clear att_sum_witha
end

for i=1:nodenum
P_with{i,k}=P0;
xe_with{i,k}=xe0;

yl_with{i,k}=lamda{i,k}*y{i,k}+(1-lamda{i,k})*H{i}*xe_with{i,k};
K_with{i,k}=(1+delta)*A*P_with{i,k}*H{i}'*inv(H{i}*P_with{i,k}*H{i}'+R1{i});
xe_with{i,k+1}=A*xe_with{i,k}+K_with{i,k}*(yl_with{i,k}-H{i}*xe_with{i,k})-attack_sum_witha{i,k};
P_with{i,k+1}=(1+delta)*(A-lamda{i,k}*K_with{i,k}*H{i})*P_with{i,k}*(A-lamda{i,k}*K_with{i,k}*H{i})'+lamda{i,k}*K_with{i,k}*R1{i}*K_with{i,k}'+(1+1/delta)*attack_sum_witha{i,k}*attack_sum_witha{i,k}'+G*Q*G';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Attack without detection
for i=1:nodenum
    for j=1:nodenum
  if L3{k}(i,j)==1 && i~=j && gama{j,k}==1
  send_xe2{i,j,k}=xe0+rou*aerfa{i,j,k};
  q2{i,j,k}=1;

  elseif without_L3{k}(i,j)==1 && i~=j && gama{j,k}==1
  send_xe2{i,j,k}=xe0;
  q2{i,j,k}=0;

  else
  send_xe2{i,j,k}=0;
  q2{i,j,k}=0;
  end
    end
end

for i=1:nodenum
    att_sum_withouta=0;
    for j=1:nodenum
        if send_xe2{i,j,k}~=0
    att_sum_withouta1=ebu*A*(1-phy{i,j,k})*(xe0-send_xe2{i,j,k});
        else
    att_sum_withouta1=0;
        end
    att_sum_withouta=att_sum_withouta+att_sum_withouta1;
    end
    attack_sum_withouta{i,k}=att_sum_withouta;
    clear att_sum_withouta
end

for i=1:nodenum
P_attack{i,k}=P0;
xe_attack{i,k}=xe0;
yl_attack{i,k}=lamda{i,k}*y{i,k}+(1-lamda{i,k})*H{i}*xe_attack{i,k};
K_attack{i,k}=(1+delta)*A*P_attack{i,k}*H{i}'*inv((1+delta)*H{i}*P_attack{i,k}*H{i}'+R1{i});
F_attack{i,k}=attack_sum_withouta{i,k};
xe_attack{i,k+1}=A*xe_attack{i,k}+K_attack{i,k}*(yl_attack{i,k}-H{i}*xe_attack{i,k})-attack_sum_withouta{i,k};
P_attack{i,k+1}=(1+delta)*(A-lamda{i,k}*K_attack{i,k}*H{i})*P_attack{i,k}*(A-lamda{i,k}*K_attack{i,k}*H{i})'+lamda{i,k}*K_attack{i,k}*R1{i}*K_attack{i,k}'+(1+1/delta)*F_attack{i,k}*F_attack{i,k}'+G*Q*G';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Attack with detection
for i=1:nodenum
    for j=1:nodenum
  if L3{k}(i,j)==1 && i~=j && lamda{j,k}==1 && gama{j,k}==1
  q{i,j,k}=1;
  z{j,k}=y{j,k}-H{j}*xe{j,k};
  send_xe{i,j,k}=xe{j,k}+rou*q{i,j,k}*aerfa{i,j,k};
  send_z{i,j,k}=z{j,k}+rou*q{i,j,k}*aerfa1{i,j,k};
  Phy{i,j,k}=1-exp(-1/2 *(send_z{i,j,k})'*inv(R1_P)*send_z{i,j,k});
     if Phy{i,j,k}>=rand1{k}
     eta{i,j,k}=0;
     else
     eta{i,j,k}=1;
     end

  elseif L3{k}(i,j)==1 && i~=j && lamda{j,k}==0 && gama{j,k}==1
  q{i,j,k}=1;
  send_xe{i,j,k}=xe{j,k}+rou*q{i,j,k}*aerfa{j,k}; 
  det_z{i,j,k}=lamda{i,k}*y{i,k}+(1-lamda{i,k})*H{i}*xe{i,k}-H{i}*send_xe{i,j,k};
  Phy{i,j,k}=1-exp(-1/2 *(det_z{i,j,k})'*inv(R1_P)*det_z{i,j,k});
     if Phy{i,j,k}>=rand1{k}
     eta{i,j,k}=0;
     else
     eta{i,j,k}=1;
     end

 elseif without_L3{k}(i,j)==1 && i~=j && lamda{j,k}==1 && gama{j,k}==1
  q{i,j,k}=0;
  z{j,k}=y{j,k}-H{j}*xe{j,k};
  send_xe{i,j,k}=xe{j,k};
  send_z {i,j,k}=z{j,k};
  Phy{i,j,k}=1-exp(-1/2 *(send_z{i,j,k})'*inv(R1_P)*send_z{i,j,k});
     if Phy{i,j,k}>=rand1{k}
     eta{i,j,k}=0;
     else
     eta{i,j,k}=1;
     end

 elseif without_L3{k}(i,j)==1 && i~=j && lamda{j,k}==0 && gama{j,k}==1
  q{i,j,k}=0;
  send_xe{i,j,k}=xe{j,k};
  det_z{i,j,k}=lamda{i,k}*y{i,k}+(1-lamda{i,k})*H{i}*xe{i,k}-H{i}*send_xe{i,j,k};
  Phy{i,j,k}=1-exp(-1/2 *(det_z{i,j,k})'*inv(R1_P)*det_z{i,j,k});
     if Phy{i,j,k}>=rand1{k}
     eta{i,j,k}=0;
     else
     eta{i,j,k}=1;
     end

  else
  q{i,j,k}=0;
  send_xe{i,j,k}=0;
  eta{i,j,k}=0;
  
  end
    end
end

for i=1:nodenum
    att_suma=0;
    for j=1:nodenum
        if  send_xe{i,j,k}~=0 & eta{i,j,k}==1
    att_suma1=ebu*A*(1-phy{i,j,k})*(xe{i,k}-send_xe{i,j,k});
        else
    att_suma1=0;
        end
    att_suma=att_suma+att_suma1;
    end
    attack_suma{i,k}=att_suma;
    clear att_suma
end


for i=1:nodenum
yl{i,k}=lamda{i,k}*y{i,k}+(1-lamda{i,k})*H{i}*xe{i,k};
K{i,k}=(1+delta)*A*P{i,k}*H{i}'*inv((1+delta)*H{i}*P{i,k}*H{i}'+R1{i});
xe{i,k+1}=A*xe{i,k}+K{i,k}*(yl{i,k}-H{i}*xe{i,k})-attack_suma{i,k};
F{i,k}=attack_suma{i,k};
P{i,k+1}=(1+delta)*(A-lamda{i,k}*K{i,k}*H{i})*P{i,k}*(A-lamda{i,k}*K{i,k}*H{i})'+lamda{i,k}*K{i,k}*R1{i}*K{i,k}'+(1+1/delta)*F{i,k}*F{i,k}'+G*Q*G';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%No attack without detection
for i=1:nodenum
  for j=1:nodenum
  if link(i,j)==1 && gama{j,k}==1
  send_xes{i,j,k}=xe_safe{j,k};
  else
  send_xes{i,j,k}=0;
  end
  end
end

for i=1:nodenum
    att_sum_safe=0;
    for j=1:nodenum 
        if send_xes{i,j,k}~=0
    att_sum_safe1=ebu*A*(xe_safe{i,k}-send_xes{i,j,k});
        else
    att_sum_safe1=0;
        end
    att_sum_safe=att_sum_safe+att_sum_safe1;
    end
    attack_sums{i,k}=att_sum_safe;
    clear att_sum_safe
end

for i=1:nodenum
yl_safe{i,k}=lamda{i,k}*y{i,k}+(1-lamda{i,k})*H{i}*xe_safe{i,k};
K_safe{i,k}=(1+delta)*A*P_safe{i,k}*H{i}'*inv((1+delta)*H{i}*P_safe{i,k}*H{i}'+R1{i});
xe_safe{i,k+1}=A*xe_safe{i,k}+K_safe{i,k}*(yl_safe{i,k}-H{i}*xe_safe{i,k})-attack_sums{i,k};
P_safe{i,k+1}=(1+delta)*(A-lamda{i,k}*K_safe{i,k}*H{i})*P_safe{i,k}*(A-lamda{i,k}*K_safe{i,k}*H{i})'+lamda{i,k}*K_safe{i,k}*R1{i}*K_safe{i,k}'+(1+1/delta)*attack_sums{i,k}*attack_sums{i,k}'+G*Q*G';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%No attack with detection
for i=1:nodenum
    for j=1:nodenum
  if link(i,j)==1 && i~=j && lamda{j,k}==1 && gama{j,k}==1
  z1{j,k}=y{j,k}-H{j}*xe_with{j,k};
  send_xe1{i,j,k}=xe_with{j,k};
  send_z1{i,j,k}=z1{j,k};
  Phy1{i,j,k}=1-exp(-1/2 *send_z1{i,j,k}'*inv(R1_P)*send_z1{i,j,k});
     if Phy1{i,j,k}>=rand1{k}
     eta1{i,j,k}=0;
     else
     eta1{i,j,k}=1;
     end

  elseif link(i,j)==1 && i~=j && lamda{j,k}==0 && gama{j,k}==1
  z1{j,k}=0;
  send_xe1{i,j,k}=xe_with{j,k}; 
  det_z1{i,j,k}=lamda{i,k}*y{i,k}-H{i}*send_xe1{i,j,k};
  Phy1{i,j,k}=1-exp(-1/2* (det_z1{i,j,k})'*inv(R1_P)*det_z1{i,j,k});
     if Phy1{i,j,k}>=rand1{k}
     eta1{i,j,k}=0;
     else
     eta1{i,j,k}=1;
     end
  
  else
  send_xe1{i,j,k}=0;
  eta1{i,j,k}=0;

  end
    end
end

for i=1:nodenum
    att_sum_witha=0;
    for j=1:nodenum
        if send_xe1{i,j,k}~=0 & eta1{i,j,k}==1
    att_sum_witha1=ebu*A*(xe_with{i,k}-send_xe1{i,j,k});
        else
    att_sum_witha1=0;
        end
        att_sum_witha=att_sum_witha+att_sum_witha1;
    end
    attack_sum_witha{i,k}=att_sum_witha;
    clear att_sum_witha
end

for i=1:nodenum
yl_with{i,k}=lamda{i,k}*y{i,k}+(1-lamda{i,k})*H{i}*xe_with{i,k};
K_with{i,k}=(1+delta)*A*P_with{i,k}*H{i}'*inv(H{i}*P_with{i,k}*H{i}'+R1{i});
xe_with{i,k+1}=A*xe_with{i,k}+K_with{i,k}*(yl_with{i,k}-H{i}*xe_with{i,k})-attack_sum_witha{i,k};
P_with{i,k+1}=(1+delta)*(A-lamda{i,k}*K_with{i,k}*H{i})*P_with{i,k}*(A-lamda{i,k}*K_with{i,k}*H{i})'+lamda{i,k}*K_with{i,k}*R1{i}*K_with{i,k}'+(1+1/delta)*attack_sum_witha{i,k}*attack_sum_witha{i,k}'+G*Q*G';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Attack without detection
for i=1:nodenum
    for j=1:nodenum
  if L3{k}(i,j)==1 && i~=j && gama{j,k}==1
  send_xe2{i,j,k}=xe_attack{j,k}+rou*aerfa{i,j,k};
  q2{i,j,k}=1;

  elseif without_L3{k}(i,j)==1 && i~=j && gama{j,k}==1
  send_xe2{i,j,k}=xe_attack{j,k};
  q2{i,j,k}=0;

  else
  send_xe2{i,j,k}=0;
  q2{i,j,k}=0;
  end
    end
end

for i=1:nodenum
    att_sum_withouta=0;
    for j=1:nodenum
        if send_xe2{i,j,k}~=0
    att_sum_withouta1=ebu*A*(1-phy{i,j,k})*(xe_attack{i,k}-send_xe2{i,j,k});
        else
            att_sum_withouta1=0;
        end
    att_sum_withouta=att_sum_withouta+att_sum_withouta1;
    end
    attack_sum_withouta{i,k}=att_sum_withouta;
    clear att_sum_withouta
end

for i=1:nodenum
yl_attack{i,k}=lamda{i,k}*y{i,k}+(1-lamda{i,k})*H{i}*xe_attack{i,k};
K_attack{i,k}=(1+delta)*A*P_attack{i,k}*H{i}'*inv((1+delta)*H{i}*P_attack{i,k}*H{i}'+R1{i});
F_attack{i,k}=attack_sum_withouta{i,k};
xe_attack{i,k+1}=A*xe_attack{i,k}+K_attack{i,k}*(yl_attack{i,k}-H{i}*xe_attack{i,k})-attack_sum_withouta{i,k};
P_attack{i,k+1}=(1+delta)*(A-lamda{i,k}*K_attack{i,k}*H{i})*P_attack{i,k}*(A-lamda{i,k}*K_attack{i,k}*H{i})'+lamda{i,k}*K_attack{i,k}*R1{i}*K_attack{i,k}'+(1+1/delta)*F_attack{i,k}*F_attack{i,k}'+G*Q*G';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%false rate
eta_1=cell2mat(eta);
eta1_1=cell2mat(eta1);

k=30;
FN(R1_P)=(length(find((L3{k}+eta_1(:,:,k))==2)))/(length(find(L3{k}==1)));
FNa(R1_P)=0;
FP(R1_P)=(length(find(without_L3{k}==1))-length(find((without_L3{k}+eta_1(:,:,k))==2)))/(length(find(without_L3{k}==1))); 
FPa(R1_P)=(length(find(link==1))-length(find((link+eta1_1(:,:,k))==2)))/(length(find(link==1))); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%RMSE
x_value=cell2mat(x);
xe_final=cell2mat(xe);

MSE_local_P=zeros(1,L+1);
MSE_with_P=zeros(1,L+1);
MSE_attack_P=zeros(1,L+1);
MSE_safe_P=zeros(1,L+1);
MSE_P=zeros(1,L+1);

MSE_local_P1=zeros(1,L+1);
MSE_with_P1=zeros(1,L+1);
MSE_attack_P1=zeros(1,L+1);
MSE_safe_P1=zeros(1,L+1);
MSE_P1=zeros(1,L+1);

for i=1:L+1
suma1=0;
suma2=0;
sumn1=0;
sumn2=0;
sumn3=0;

xea_final=xe(:,i);
xea_safe_final=xe_safe(:,i);
xea_with_final=xe_with(:,i);
xea_attack_final=xe_attack(:,i);

for j=1:nodenum
suma1=suma1+(xea_final{j,:}-x_value(:,i)).*(xea_final{j,:}-x_value(:,i));
suma2=suma2+(xea_safe_final{j,:}-x_value(:,i)).*(xea_safe_final{j,:}-x_value(:,i));
sumn1=sumn1+(xea_with_final{j,:}-x_value(:,i)).*(xea_with_final{j,:}-x_value(:,i));
sumn2=sumn2+(xea_attack_final{j,:}-x_value(:,i)).*(xea_attack_final{j,:}-x_value(:,i));
end

suma_1=suma1(1,:)+suma1(3,:); 
suma_2=suma2(1,:)+suma2(3,:); 
sum1=sumn1(1,:)+sumn1(3,:); 
sum2=sumn2(1,:)+sumn2(3,:); 

suma_1v=suma1(2,:)+suma1(4,:); 
suma_2v=suma2(2,:)+suma2(4,:);
sum1v=sumn1(2,:)+sumn1(4,:); 
sum2v=sumn2(2,:)+sumn2(4,:); 

MSE1_pos(1,i,Mont)=suma_1;
MSE2_pos(1,i,Mont)=suma_2;
MSE3_pos(1,i,Mont)=sum1;
MSEa_pos(1,i,Mont)=sum2;

MSE1_vel(1,i,Mont)=suma_1v;
MSE2_vel(1,i,Mont)=suma_2v;
MSE3_vel(1,i,Mont)=sum1v;
MSEa_vel(1,i,Mont)=sum2v;
end
end

for s=1:L+1
sumb=0;
sumc=0;
sumd=0;
sume=0;

    for Mont=1:MC
sumb=sumb+MSE1_pos(1,s,Mont);
sumc=sumc+MSE2_pos(1,s,Mont);
sumd=sumd+MSE3_pos(1,s,Mont);
sume=sume+MSEa_pos(1,s,Mont);

sumb1=sumb+MSE1_vel(1,s,Mont);
sumc1=sumc+MSE2_vel(1,s,Mont);
sumd1=sumd+MSE3_vel(1,s,Mont);
sume1=sume+MSEa_vel(1,s,Mont);
    end
    MSE_P(s)=sumb/MC/nodenum;
    MSE_safe_P(s)=sumc/MC/nodenum;
    MSE_with_P(s)=sumd/MC/nodenum;
    MSE_attack_P(s)=sume/MC/nodenum; 

    MSE_P1(s)=sumb1/MC/nodenum;
    MSE_safe_P1(s)=sumc1/MC/nodenum;
    MSE_with_P1(s)=sumd1/MC/nodenum;
    MSE_attack_P1(s)=sume1/MC/nodenum;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_P_F=zeros(1,L+1);
sum_P_F1=zeros(1,L+1);
sum_P_F_without=zeros(1,L+1);
sum_P_F_without1=zeros(1,L+1);
sum_P_F_with=zeros(1,L+1);
sum_P_F_with1=zeros(1,L+1);
sum_P_F_without_A=zeros(1,L+1);
sum_P_F_without_A1=zeros(1,L+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:L+1
    P1{i}=P(:,i);
    sum_P=0;
    sum_P1=0;
    for j=1:nodenum
        sum_P=sum_P+P1{i}{j,1}(1,1)+P1{i}{j,1}(3,3);
        sum_P1=sum_P1+P1{i}{j,1}(2,2)+P1{i}{j,1}(4,4);
    end
    sum_P_F(1,i)=sum_P/(nodenum);
    sum_P_F1(1,i)=sum_P1/(nodenum);
    clear sum_P
    clear sum_P1
end

for i=1:L+1
    P1_without{i}=P_safe(:,i);
    sum_P=0;
    sum_P1=0;
    for j=1:nodenum
        sum_P=sum_P+P1_without{i}{j,1}(1,1)+P1_without{i}{j,1}(3,3);
        sum_P1=sum_P1+P1_without{i}{j,1}(2,2)+P1_without{i}{j,1}(4,4);
    end
    sum_P_F_without(1,i)=sum_P/(nodenum);
    sum_P_F_without1(1,i)=sum_P1/(nodenum);
    clear sum_P
    clear sum_P1
end

for i=1:L+1
    P1_with{i}=P_with(:,i);
    sum_P=0;
    sum_P1=0;
    for j=1:nodenum
        sum_P=sum_P+P1_with{i}{j,1}(1,1)+P1_with{i}{j,1}(3,3);
        sum_P1=sum_P1+P1_with{i}{j,1}(2,2)+P1_with{i}{j,1}(4,4);
    end
    sum_P_F_with(1,i)=sum_P/(nodenum);
    sum_P_F_with1(1,i)=sum_P1/(nodenum);
    clear sum_P
    clear sum_P1
end

for i=1:L+1
    P1_without_A{i}=P_attack(:,i);
    sum_P=0;
    sum_P1=0;
    for j=1:nodenum
        sum_P=sum_P+P1_without_A{i}{j,1}(1,1)+P1_without_A{i}{j,1}(3,3);
        sum_P1=sum_P1+P1_without_A{i}{j,1}(2,2)+P1_without_A{i}{j,1}(4,4);
    end
    sum_P_F_without_A(1,i)=sum_P/(nodenum);
    sum_P_F_without_A1(1,i)=sum_P1/(nodenum);
    clear sum_P
    clear sum_P1
end

jz=zeros(1,nodenum);
jz1=zeros(1,nodenum);
for i=0:1:nodenum-1
    jz(i+1)=4*i+1;
    jz1(i+1)=4*i+2;
end
xe1_final=xe_final(jz,:);
xe2_final=xe_final(jz1,:);

figure(1);
hold on 
p3=plot(x_value(1,:),x_value(3,:),'k','linewidth',2);
legend([p1,p3],{'Sensor','Target trajectory'},'FontName','Arial');
xlabel('Coordinate X')
ylabel('Coordinate Y')
xlim([0, 100])
set(gca,'FontName','Arial');
 
figure(2)
plot(1:L+1,sum_P_F(1,:),'c-*',1:L+1,sum_P_F_without(1,:),'b-.',1:L+1,sum_P_F_with(1,:),'r-',1:L+1,sum_P_F_without_A(1,:),'k-x','linewidth',1)
legend({'1. Attack with detection','2. No attack without detection','3. No attack with detection','4. Attack without detection'},'FontName','Arial');
xlabel('Time')
ylabel('J_k^{Pos}')
set(gca,'FontName','Arial');

figure(3)
plot(1:L+1,sum_P_F1(1,:),'c-*',1:L+1,sum_P_F_without1(1,:),'b-.',1:L+1,sum_P_F_with1(1,:),'r-',1:L+1,sum_P_F_without_A1(1,:),'k-x','linewidth',1)
legend({'1. Attack with detection','2. No attack without detection','3. No attack with detection','4. Attack without detection'},'FontName','Arial');
xlabel('Time')
ylabel('J_k^{Vel}')
set(gca,'FontName','Arial');


figure(4)
plot(1:L+1,x_value(1,:),'k',1:L+1,xe1_final(1:nodenum,:),'linewidth',1)
legend({'Target','Sensor 1-20'},'FontName','Arial');
xlabel('Time')
ylabel('The positions on the X-axis')
xlim([0, 100])
set(gca,'FontName','Arial');
 
figure(5)
plot(1:L+1,x_value(2,:),'k',1:L+1,xe2_final(1:nodenum,:),'linewidth',1)
legend({'Target','Sensor 1-20'},'FontName','Arial');
xlabel('Time')
ylabel('The velocities on the X-axis')
xlim([0, 100])
set(gca,'FontName','Arial');


figure(6)
plot(1:L+1,MSE_P(1,:),'c-*',1:L+1,MSE_safe_P(1,:),'b-.',1:L+1,MSE_with_P(1,:),'r-',1:L+1,MSE_attack_P(1,:),'k-x','linewidth',1)
legend({'1. Attack with detection','2. No attack without detection','3. No attack with detection','4. Attack without detection'},'FontName','Arial');
xlabel('Time')
ylabel('MSE_{k}^{pos}')
xlim([0, 100])
ylim([0, 3])
set(gca,'FontName','Arial');

figure(7)
plot(1:L+1,MSE_P1(1,:),'c-*',1:L+1,MSE_safe_P1(1,:),'b-.',1:L+1,MSE_with_P1(1,:),'r-',1:L+1,MSE_attack_P1(1,:),'k-x','linewidth',1)
legend({'1. Attack with detection','2. No attack without detection','3. No attack with detection','4. Attack without detection'},'FontName','Arial');
xlabel('Time')
ylabel('MSE_{k}^{vel}')
xlim([0, 100])
ylim([0, 3])
set(gca,'FontName','Arial');
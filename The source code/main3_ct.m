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

%Initialization state, estimated value and variance
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
MC=10;
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
for R1_P=1:LP
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
%Dos attak
dosProb=0.2;
L1{k}=binornd(1,dosProb,nodenum).*link;
L2{k}=triu(L1{k})+triu(L1{k})';
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
L3{k}=binornd(1,fdiProb,nodenum).*link; 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
end
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fale rate
eta_1=cell2mat(eta);
k=30;
FN(R1_P)=(length(find((L3{k}+eta_1(:,:,k))==2)))/(length(find(L3{k}==1)));
FP(R1_P)=(length(find(without_L3{k}==1))-length(find((without_L3{k}+eta_1(:,:,k))==2)))/(length(find(without_L3{k}==1))); 
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

for j=1:nodenum
suma1=suma1+(xea_final{j,:}-x_value(:,i)).*(xea_final{j,:}-x_value(:,i));
end

suma_1=suma1(1,:)+suma1(3,:); 
suma_1v=suma1(2,:)+suma1(4,:);

MSE1_pos(1,i,Mont)=suma_1;
MSE1_vel(1,i,Mont)=suma_1v;
end
end

for s=1:L+1
sumb=0;
sumc=0;
sumd=0;
sume=0;
sumf=0;

for Mont=1:MC
sumb=sumb+MSE1_pos(1,s,Mont);
sumc=sumc+MSE2_pos(1,s,Mont);
sumd=sumd+MSE3_pos(1,s,Mont);
sume=sume+MSEa_pos(1,s,Mont);
sumf=sumf+MSEb_pos(1,s,Mont);

sumb1=sumb+MSE1_vel(1,s,Mont);
sumc1=sumc+MSE2_vel(1,s,Mont);
sumd1=sumd+MSE3_vel(1,s,Mont);
sume1=sume+MSEa_vel(1,s,Mont);
sumf1=sumf+MSEb_vel(1,s,Mont);
end
    MSE_P(s)=sumb/MC/nodenum;
    MSE_safe_P(s)=sumc/MC/nodenum;
    MSE_with_P(s)=sumd/MC/nodenum;
    MSE_attack_P(s)=sume/MC/nodenum;
    MSE_local_P(s)=sumf/MC/nodenum; 

    MSE_P1(s)=sumb1/MC/nodenum;
    MSE_safe_P1(s)=sumc1/MC/nodenum;
    MSE_with_P1(s)=sumd1/MC/nodenum;
    MSE_attack_P1(s)=sume1/MC/nodenum;
    MSE_local_P1(s)=sumf1/MC/nodenum;
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
sum_P1_local=zeros(1,L+1);
sum_P1_local1=zeros(1,L+1);

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

figure(1);
hold on 
p3=plot(x_value(1,:),x_value(3,:),'k','linewidth',2);
legend([p1,p3],{'Sensor','Target trajectory'},'FontName','Arial');
xlabel('Coordinate X')
ylabel('Coordinate Y')
set(gca,'FontName','Arial');

x=1:LP;
y=polyfit(x,FP(1,:),2);
y1=polyfit(x,FN(1,:),2);

ya=polyfit(x,FPa(1,:),2);
ya1=polyfit(x,FNa(1,:),2);
z=polyval(y,x);
z1=polyval(y1,x);

za=polyval(ya,x);
za1=polyval(ya1,x);

figure(2)
plot(x,z,'r-*',x,z1,'b-+')
legend({'FP','FN'},'FontName','Arial');
xlabel('Parameter Z')
ylabel('False rate')
set(gca,'FontName','Arial');

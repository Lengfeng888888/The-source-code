clear;
numofnodes = 20;

figure(1);
clf;
xticks(0:20:100);
yticks(0:20:100);
R = 30; % link range
Ra = 20; % link range
sen1_x = readmatrix('network_attack.xlsx','Sheet','Sheet1','Range','A3:A12');
sen1_y = readmatrix('network_attack.xlsx','Sheet','Sheet1','Range','B3:B12');
sen2_x = readmatrix('network_attack.xlsx','Sheet','Sheet1','Range','C3:C12');
sen2_y = readmatrix('network_attack.xlsx','Sheet','Sheet1','Range','D3:D12');
link=zeros(numofnodes);
netx = [sen1_x;sen2_x];
nety = [sen1_y;sen2_y];
net=[netx nety];


scatter(netx,nety,'filled');
p1 = plot(netx,nety,'o','MarkerSize',8,'MarkerEdgeColor','r','MarkerFaceColor','r');
hold on;

for i = 1 : numofnodes
    for j = 1 : numofnodes
        distance = sqrt((netx(i) - netx(j))^2 + (nety(i) - nety(j))^2);
        if (distance <= R && i~=j && i<j)
            link(i,j) = 1;
            draw_arrow([netx(i),nety(i)],[netx(j),nety(j)]);
        elseif  (distance <= Ra && i~=j && i>j)
             link(i,j) = 1;
            draw_arrow([netx(i),nety(i)],[netx(j),nety(j)]);
        else 
            link(i,j) = 0;
        end
    end
end

text(net(:,1),net(:,2),arrayfun(@(net)['  ' num2str(net)],1:numofnodes,'UniformOutput',0))

hold on;
figure(1);
xlim([0 105])
ylim([0 105])

legend(p1,{'Sensor'},'FontName','Arial');
xlabel('Coordinate x [m]','FontName','Arial')
ylabel('Coordinate y [m]','FontName','Arial')
function draw_arrow(start_point, end_point)
K = 0.06;  
theta = pi / 8; 
A1 = [cos(theta), -sin(theta);
     sin(theta), cos(theta)]; 
theta = -theta;
A2 = [cos(theta), -sin(theta);
     sin(theta), cos(theta)]; 
arrow = start_point' - end_point';
arrow_1 = A1 * arrow;
arrow_2 = A2 * arrow;
arrow_1 = K * arrow_1 + end_point';
arrow_2 = K * arrow_2 + end_point';

plot([start_point(1), end_point(1)], [start_point(2), end_point(2)], 'b','LineWidth',1);
plot([arrow_1(1), end_point(1)], [arrow_1(2), end_point(2)], 'b','LineWidth',1);
plot([arrow_2(1), end_point(1)], [arrow_2(2), end_point(2)], 'b','LineWidth',1);
end








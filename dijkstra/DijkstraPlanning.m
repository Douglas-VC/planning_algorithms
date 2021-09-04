close all;
clear;
tic

%% Parâmetros de Entrada

% metric = 'slope';
metric = 'distance';
% metric = 'height';
% metric = 'combined';

adjustFlatTerrain = 0;

dec = 3;
dec_terreno = 3;

addpath('./Libraries/graphutils');

ptCloud = pcread('./PointClouds/campinhoremaster3.pcd');
% ptCloud = pcread('./PointClouds/campinhoremaster3reduzido.pcd');
% ptCloud = pcread('./PointClouds/map_lego_reduzido_filtrado.ply');
% ptCloud = pcread('./PointClouds/flatTerrain.ply');

figure (2)
pcshow(ptCloud);
view(15,25)
grid off
set(gca,'visible','off')
set(gcf, 'InvertHardCopy', 'off');
set(gcf,'Color',[0 0 0]);
set(gca,'FontSize',30)
xlabel('x')
ylabel('y')
zlabel('z')

p_start = [13.8 -4.1];
p_obj = [-15 1.25];

% p_start = [4 4];
% p_obj = [3 -4];
% 
% p_start = [3 -4];
% p_obj = [4 4];


%% Processamento dos dados da point cloud e geração de função h

xyz = ptCloud.Location;
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

[Uxy,~, idx]  = unique([x y], 'rows');
zmax = accumarray(idx, z, [], @max);
XYZ = [Uxy zmax];

X = round(XYZ(:,1),4);
Y = round(XYZ(:,2),4);
Z = round(XYZ(:,3),4);

step = 0.1;
x_grid = min(X) : step : max(X);
y_grid = min(Y) : step : max(Y);

[X2,Y2] = meshgrid(x_grid,y_grid);
X2 = round(X2,4);
Y2 = round(Y2,4);
Z2 = griddata(double(X),double(Y),double(Z),double(X2),double(Y2),'cubic');
Z2(isnan(Z2)) = 0;

[b,a] = butter(3,0.80);
for k = 1:1:size(Z2,1)
    Z2(k,:) = filtfilt(b,a,Z2(k,:));
end
for k = 1:1:size(Z2,2)
    Z2(:,k) = filtfilt(b,a,Z2(:,k));
end

X3 = X2(1:dec:end,1:dec:end);
Y3 = Y2(1:dec:end,1:dec:end);
Z3 = Z2(1:dec:end,1:dec:end);

X6 = X2(1:dec_terreno:end,1:dec_terreno:end);
Y6 = Y2(1:dec_terreno:end,1:dec_terreno:end);
Z6 = Z2(1:dec_terreno:end,1:dec_terreno:end);

global function_h;
function_h = griddedInterpolant(X6',Y6',Z6');
function_h.Method = 'spline';

%% Plot do terreno

x4 = min(X2(:)):0.05:max(X2(:));
y4 = min(Y2(:)):0.05:max(Y2(:));
[X4,Y4] = meshgrid(x4,y4);
Z4 = function_h(X4',Y4');
Z4 = Z4';

% Plot the interpolation result X4
% figure(101)
% handle = surf(X4,Y4,Z4);
% handle.LineStyle = 'none';
% axis equal
% grid on
% xlabel('$x$','Interpreter','latex','FontSize',15)
% ylabel('$y$','Interpreter','latex','FontSize',15)
% zlabel('$z$','Interpreter','latex','FontSize',15)
% title('X4','FontSize',15)
% 
% % Plot the interpolation result X3
% figure(102)
% handle = surf(X3,Y3,Z3);
% handle.LineStyle = 'none';
% axis equal
% grid on
% xlabel('$x$','Interpreter','latex','FontSize',15)
% ylabel('$y$','Interpreter','latex','FontSize',15)
% zlabel('$z$','Interpreter','latex','FontSize',15)
% title('X3','FontSize',15)

%% Organização das coordenadas dos vértices presentes no grafo

N = (size(X3(1,:),2)*size(X3(:,1),1));
coordinates = zeros(N,2);
C = zeros(N,N,'single');

index = 1;
for j = 1:size(X3(:,1),1)
    for i = 1:size(X3(1,:),2) 
        coordinates(index,1) = X3(1,i);
        coordinates(index,2) = Y3(end - j + 1,1);
        index = index + 1;
    end
end

%% Criação da matriz de adjacência A que representa as conexões do grafo

x_size = size(X3(1,:),2);
y_size = size(Y3(:,1),1);
A = C*0;

for i = 1:N
    for j = i:N
        if (j == i+1 && mod(i,x_size) ~= 0) || j == i+x_size || (j == i+x_size-1 && mod(j,x_size) ~= 0) || (j == i+x_size+1 && mod(i,x_size) ~= 0)
%         if (j == i+1 && mod(i,x_size) ~= 0) || j == i+x_size

            A(i,j) = 1;
            A(j,i) = 1;
            
        end
    end
end

% figure (2)
% hold on
% % scatter(coordinates(:,1),coordinates(:,2),'filled', 'k');
% scatter3(coordinates(:,1),coordinates(:,2),function_h(coordinates(:,1),coordinates(:,2)),'filled', 'w');
% for i = 1:N
%     for j = i:N
%         if A(i,j) == 1
% %             plot([coordinates(i,1) coordinates(j,1)],[coordinates(i,2) coordinates(j,2)],'k');
%             plot3([coordinates(i,1) coordinates(j,1)],[coordinates(i,2) coordinates(j,2)],[function_h(coordinates(i,1),coordinates(i,2))+0.1,function_h(coordinates(j,1),coordinates(j,2))+0.1], 'w');
%         end
%     end
% end
% % hold off
% set(gca,'visible','off')

%% Cálculo da matriz de custos C a partir da métrica escolhida

X5 = X3;
Y5 = Y3;
Z5 = X5*inf;

if strcmp(metric,'slope')
    delta = 1e-3;
    for i = 1:1:length(X5(1,:))
        for j = 1:1:length(X5(:,1))

            p = [X5(j,i); Y5(j,i)];        

            Pa = [p(1)       p(2)+delta function_h([p(1), p(2)+delta])];
            Pb = [p(1)+delta p(2)-delta function_h([p(1)+delta, p(2)-delta])];
            Pc = [p(1)-delta p(2)-delta function_h([p(1)-delta, p(2)-delta])];

            normal = cross(Pa-Pb, Pa-Pc);
            angle = acosd(normal * [0 0 -1]'/norm(normal));

            Z5(j,i) = angle;
        end
    end
    function_c = griddedInterpolant(X5',Y5',Z5');
    function_c.Method = 'spline';
    
    for i = 1:N
        for j = 1:N
            if A(i,j) == 1
                p1 = coordinates(i,:);
                p2 = coordinates(j,:);
                cost = abs((function_c(p1)+function_c(p2)+function_c((p1+p2)/2))/3 + max([function_c(p1),function_c(p2),function_c((p1+p2)/2)]));
                C(i,j) = cost;
                C(j,i) = cost;
            end
        end
    end
    
elseif strcmp(metric,'height')
    
    for i = 1:1:length(X5(1,:))
        for j = 1:1:length(X5(:,1))
            Z5(j,i) = f_cost([X5(j,i); Y5(j,i)]);
        end
    end
    
    function_c = griddedInterpolant(X5',Y5',Z5');
    function_c.Method = 'spline';
    
    for i = 1:N
        for j = 1:N
            if A(i,j) == 1
                p1 = coordinates(i,:);
                p2 = coordinates(j,:);
                dist = sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2);
                f = (function_c(p1) + function_c((p1+p2)/2) + function_c(p2))/3;
                cost = f*dist;
                C(i,j) = cost;
                C(j,i) = cost;
            end
        end
    end
    
elseif strcmp(metric,'combined')
    
    for i = 1:1:length(X5(1,:))
        for j = 1:1:length(X5(:,1))
            Z5(j,i) = f_cost([X5(j,i); Y5(j,i)]);
        end
    end
    
    function_c_height = griddedInterpolant(X5',Y5',Z5');
    function_c_height.Method = 'spline';
    
    delta = 1e-3;
    for i = 1:1:length(X5(1,:))
        for j = 1:1:length(X5(:,1))

            p = [X5(j,i); Y5(j,i)];        

            Pa = [p(1)       p(2)+delta function_h([p(1), p(2)+delta])];
            Pb = [p(1)+delta p(2)-delta function_h([p(1)+delta, p(2)-delta])];
            Pc = [p(1)-delta p(2)-delta function_h([p(1)-delta, p(2)-delta])];

            normal = cross(Pa-Pb, Pa-Pc);
            angle = acosd(normal * [0 0 -1]'/norm(normal));

            Z5(j,i) = angle;
        end
    end
    function_c_inc = griddedInterpolant(X5',Y5',Z5');
    function_c_inc.Method = 'spline';
    
    k1 = 6.5; %Peso para métrica de distância
    k2 = 20;   %Peso para métrica de altura
    k3 = 1/20; %Peso para métrica de inclinação
    
    for i = 1:N
        for j = 1:N
            if A(i,j) == 1
                p1 = coordinates(i,:);
                p2 = coordinates(j,:);
                inc = abs((function_c_inc(p1)+function_c_inc(p2)+function_c_inc((p1+p2)/2))/3 + max([function_c_inc(p1),function_c_inc(p2),function_c_inc((p1+p2)/2)]));
                dist = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));
                f = (function_c_height(p1) + function_c_height((p1+p2)/2) + function_c_height(p2))/3;
                cost = abs(k1*dist + k2*f*dist + k3*inc);
                C(i,j) = cost;
                C(j,i) = cost;
            end
        end
    end

elseif strcmp(metric,'distance')
    
    for i = 1:N
        for j = 1:N
            if A(i,j) == 1
                p1 = coordinates(i,:);
                p2 = coordinates(j,:);
                cost = sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2);
%                 cost = sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2);
                C(i,j) = cost;
                C(j,i) = cost;
            end
        end
    end
    
end

if adjustFlatTerrain && (strcmp(metric,'slope') || strcmp(metric,'height'))
    C = AdjustCostforFlatTerrain(x_size,y_size,dec,function_h,coordinates,C);
end
%% Determinação dos vértices inicial e final a partir dos pontos inicial e final desejados

if ismembertol(p_start,coordinates,'ByRows',0.0001)
    [~,ID_start] = ismembertol(p_start,coordinates,'ByRows',0.0001);
else
    dist_x_min_y_min = inf;
    dist_x_min_y_max = inf;
    dist_x_max_y_min = inf;
    dist_x_max_y_max = inf;
    for i = 1:length(coordinates)
        if coordinates(i,1) < p_start(1) && coordinates(i,2) < p_start(2)
            if norm(coordinates(i,:)-p_start) < dist_x_min_y_min
                start_points(1,:) = coordinates(i,:);
                dist_x_min_y_min = norm(coordinates(i,:)-p_start);
            end
        end
        if coordinates(i,1) < p_start(1) && coordinates(i,2) > p_start(2)
            if norm(coordinates(i,:)-p_start) < dist_x_min_y_max
                start_points(2,:) = coordinates(i,:);
                dist_x_min_y_max = norm(coordinates(i,:)-p_start);
            end
        end
        if coordinates(i,1) > p_start(1) && coordinates(i,2) < p_start(2)
            if norm(coordinates(i,:)-p_start) < dist_x_max_y_min
                start_points(3,:) = coordinates(i,:);
                dist_x_max_y_min = norm(coordinates(i,:)-p_start);
            end
        end
        if coordinates(i,1) > p_start(1) && coordinates(i,2) > p_start(2)
            if norm(coordinates(i,:)-p_start) < dist_x_max_y_max
                start_points(4,:) = coordinates(i,:);
                dist_x_max_y_max = norm(coordinates(i,:)-p_start);
            end
        end
    end  
    distances_from_objective = sqrt(sum(bsxfun(@minus, start_points, p_obj).^2,2));
    [~,I] = min(distances_from_objective);
    ID_start = find(ismember(coordinates, start_points(I,:),'rows'));
end

if ismembertol(p_obj,coordinates,'ByRows',0.0001)
    [~,ID_obj] = ismembertol(p_obj,coordinates,'ByRows',0.0001);
else
    dist_x_min_y_min = inf;
    dist_x_min_y_max = inf;
    dist_x_max_y_min = inf;
    dist_x_max_y_max = inf;
    for i = 1:length(coordinates)
        if coordinates(i,1) < p_obj(1) && coordinates(i,2) < p_obj(2)
            if norm(coordinates(i,:)-p_obj) < dist_x_min_y_min
                obj_points(1,:) = coordinates(i,:);
                dist_x_min_y_min = norm(coordinates(i,:)-p_obj);
            end
        end
        if coordinates(i,1) < p_obj(1) && coordinates(i,2) > p_obj(2)
            if norm(coordinates(i,:)-p_obj) < dist_x_min_y_max
                obj_points(2,:) = coordinates(i,:);
                dist_x_min_y_max = norm(coordinates(i,:)-p_obj);
            end
        end
        if coordinates(i,1) > p_obj(1) && coordinates(i,2) < p_obj(2)
            if norm(coordinates(i,:)-p_obj) < dist_x_max_y_min
                obj_points(3,:) = coordinates(i,:);
                dist_x_max_y_min = norm(coordinates(i,:)-p_obj);
            end
        end
        if coordinates(i,1) > p_obj(1) && coordinates(i,2) > p_obj(2)
            if norm(coordinates(i,:)-p_obj) < dist_x_max_y_max
                obj_points(4,:) = coordinates(i,:);
                dist_x_max_y_max = norm(coordinates(i,:)-p_obj);
            end
        end
    end
    distances_from_start = sqrt(sum(bsxfun(@minus, obj_points, p_start).^2,2));
    [~,I] = min(distances_from_start);
    ID_obj = find(ismember(coordinates, obj_points(I,:),'rows'));
end

%% Algoritmo Dijkstra que recebe matriz A de vértices do grafo, matriz C com custos, vértice inicial e vértice final

[costs,paths] = dijkstra(A,C,ID_start,ID_obj,0);

%% Processamento da trajetória gerada

pts = [];
for index = 1:size(paths(1,:),2)
    pts = [pts, [coordinates(paths(index),1) coordinates(paths(index),2)]'];
end

pts2 = pts(:,1);
n_insert = 10;
alpha = linspace(1/n_insert,1,n_insert);
% alpha = [1; 1]*alpha;
for k = 1:1:length(pts(1,:))-1
    new_pts = pts(:,k)*(1-alpha) + pts(:,k+1)*(alpha);
    pts2 = [pts2, new_pts];
end
pts2 = average_filter(pts2,10,1);
for k = 1:1:length(pts2(1,:))
    pts2(3,k) = function_h(pts2(1,k),pts2(2,k));
end

%% Adição de pontos inicial e final à trajetória
if ~ismembertol(p_start,coordinates,'ByRows',0.0001)
    pts = [p_start' pts];
    p1 = p_start;
    p2 = pts(:,2)';
    if strcmp(metric,'slope')
        cost = abs((function_c(p1)+function_c(p2)+function_c((p1+p2)/2))/3 + max([function_c(p1),function_c(p2),function_c((p1+p2)/2)]));
    elseif strcmp(metric,'heigth')
        dist = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));
        f = (function_c(p1) + function_c((p1+p2)/2) + function_c(p2))/3;
        cost = abs(f*dist);
    elseif strcmp(metric,'combined')
        inc = abs((function_c_inc(p1)+function_c_inc(p2)+function_c_inc((p1+p2)/2))/3 + max([function_c_inc(p1),function_c_inc(p2),function_c_inc((p1+p2)/2)]));
        dist = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));
        f = (function_c_height(p1) + function_c_height((p1+p2)/2) + function_c_height(p2))/3;
        cost = abs(k1*dist + k2*f*dist + k3*inc);
    elseif strcmp(metric,'distance')
        cost = sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2);
    end
    costs = costs + cost;
end
if ~ismembertol(p_obj,coordinates,'ByRows',0.0001)
    pts = [pts p_obj'];
    p1 = p_obj;
    p2 = pts(:,end-1)';
    if strcmp(metric,'slope')
        cost = abs((function_c(p1)+function_c(p2)+function_c((p1+p2)/2))/3 + max([function_c(p1),function_c(p2),function_c((p1+p2)/2)]));
    elseif strcmp(metric,'heigth')
        dist = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));
        f = (function_c(p1) + function_c((p1+p2)/2) + function_c(p2))/3;
        cost = abs(f*dist);
    elseif strcmp(metric,'combined')
        inc = abs((function_c_inc(p1)+function_c_inc(p2)+function_c_inc((p1+p2)/2))/3 + max([function_c_inc(p1),function_c_inc(p2),function_c_inc((p1+p2)/2)]));
        dist = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));
        f = (function_c_height(p1) + function_c_height((p1+p2)/2) + function_c_height(p2))/3;
        cost = abs(k1*dist + k2*f*dist + k3*inc);
    elseif strcmp(metric,'distance')
        cost = sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2);
    end
    costs = costs + cost;
end

%% Plot do terreno e trajetória gerada

figure(1)
handle = surf(X4,Y4,Z4);
handle.LineStyle = 'none';
hold on
% 
% % scatter(coordinates(:,1),coordinates(:,2),'filled', 'k');
% scatter3(coordinates(:,1),coordinates(:,2),function_h(coordinates(:,1),coordinates(:,2)),'filled', 'k');
% for i = 1:N
%     for j = i:N
%         if A(i,j) == 1
% %             plot([coordinates(i,1) coordinates(j,1)],[coordinates(i,2) coordinates(j,2)],'k');
%             plot3([coordinates(i,1) coordinates(j,1)],[coordinates(i,2) coordinates(j,2)],[function_h(coordinates(i,1),coordinates(i,2))+0.1,function_h(coordinates(j,1),coordinates(j,2))+0.1], 'k');
%         end
%     end
% end

% plot3(pts2(1,:),pts2(2,:),pts2(3,:)+0.01,'r','LineWidth',3)
plot3(pts(1,:),pts(2,:),function_h(pts(1,:),pts(2,:))+0.1,'r','LineWidth',3);
scatter3(pts(1,:),pts(2,:),function_h(pts(1,:),pts(2,:))+0.1,'o');

hold off
axis equal

fprintf("\nDistância entre vértices do grafo: %d cm\n", dec*10);
fprintf("Custo normalizado: %.3f\n", costs/size(paths,2));
fprintf("Custo total: %.3f\n", costs);
fprintf("Tempo de Processamento: %.3f segundos\n", toc);
%%
name = './logs/';
timestamp = now;
timestamp = datetime(timestamp,'ConvertFrom','datenum','Format','d-MMM-y HH-mm-ss');
timestamp = sprintf(' %s',timestamp);
name = strcat(name,sprintf('%s',metric));
name = strcat(name,timestamp);
name = strcat(name,'.mat');
save(name,'p_start','p_obj','X4','Y4','Z4','function_h','pts','pts2');

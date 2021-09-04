function C = AdjustCostforFlatTerrain(x_size,y_size,dec,function_h,coordinates,C)

graphFlatTerrain = zeros(y_size, x_size);
line_index = 1;
column_index = 1;
for ID_number = 1:length(coordinates)
    cima = 1;
    baixo = 1;
    esquerda = 1;
    direita = 1;

    if line_index > 1 && C(ID_number, ID_number - x_size) > 0.001
        cima = 0;
    end
    if line_index < y_size && C(ID_number, ID_number + x_size) > 0.001
        baixo = 0;
    end
    if column_index > 1 && C(ID_number, ID_number - 1) > 0.001
        esquerda = 0;
    end
    if column_index < x_size && C(ID_number, ID_number + 1) > 0.001
        direita = 0;
    end

    if cima && baixo && esquerda && direita
        graphFlatTerrain(line_index, column_index) = 1;
    end
    if column_index == x_size
        line_index = line_index + 1;
        column_index = 1;
    else
        column_index = column_index + 1; 
    end      

end

flatTerrainRegions = regionprops(logical(graphFlatTerrain),'boundingbox');
for i = 1:length(flatTerrainRegions)
    if ((flatTerrainRegions(i).BoundingBox(3)-1)*dec/10)*((flatTerrainRegions(i).BoundingBox(4)-1)*dec/10) > 2
        
        line_start = ceil(flatTerrainRegions(i).BoundingBox(2));
        column_start = ceil(flatTerrainRegions(i).BoundingBox(1));
        line_finish = line_start + flatTerrainRegions(i).BoundingBox(4) - 1;
        column_finish = column_start + flatTerrainRegions(i).BoundingBox(3) - 1;
        
        ID_number_array = [];
        
        for j = line_start:line_finish
            for k = column_start:column_finish
                if graphFlatTerrain(j,k) == 1
                    ID_number_array = [ID_number_array (j - 1)*x_size + k];
                end
            end
        end
        
        for j = 1:length(ID_number_array)
            
            p1 = coordinates(ID_number_array(j),:);
            
            if ID_number_array(j) > x_size
                p2 = coordinates(ID_number_array(j) - x_size,:);
                C(ID_number_array(j), ID_number_array(j) - x_size) = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));
                C(ID_number_array(j)- x_size, ID_number_array(j)) = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));

            end
            
            if ID_number_array(j) < (length(coordinates)+ 1 - x_size)
                p2 = coordinates(ID_number_array(j) + x_size,:);
                C(ID_number_array(j), ID_number_array(j) + x_size) = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));
                C(ID_number_array(j)+ x_size, ID_number_array(j)) = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));

            end
            
            if mod(ID_number_array(j), x_size) ~= 0
                p2 = coordinates(ID_number_array(j) + 1,:);
                C(ID_number_array(j), ID_number_array(j) + 1) = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));
                C(ID_number_array(j) + 1, ID_number_array(j)) = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));

            end
            
            if mod(ID_number_array(j), x_size+1) ~= 0 && ID_number_array(j) ~= 1
                p2 = coordinates(ID_number_array(j) - 1,:);
                C(ID_number_array(j), ID_number_array(j) - 1) = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));
                C(ID_number_array(j) - 1, ID_number_array(j)) = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));

            end
            
            if ID_number_array(j) > x_size && mod(ID_number_array(j), x_size) ~= 0
                p2 = coordinates(ID_number_array(j) - x_size + 1,:);
                C(ID_number_array(j), ID_number_array(j) - x_size + 1) = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));
                C(ID_number_array(j) - x_size + 1, ID_number_array(j)) = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));

            end
            
            if ID_number_array(j) > x_size && mod(ID_number_array(j), x_size+1) ~= 0 && ID_number_array(j) ~= 1
                p2 = coordinates(ID_number_array(j) - x_size - 1,:);
                C(ID_number_array(j), ID_number_array(j) - x_size - 1) = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));
                C(ID_number_array(j) - x_size - 1, ID_number_array(j)) = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));

            end
            
            if ID_number_array(j) < (length(coordinates)+ 1 - x_size)  && mod(ID_number_array(j), x_size) ~= 0
                p2 = coordinates(ID_number_array(j) + x_size + 1,:);
                C(ID_number_array(j), ID_number_array(j) + x_size + 1) = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));
                C(ID_number_array(j) + x_size + 1, ID_number_array(j)) = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));

            end
            
            if ID_number_array(j) < (length(coordinates)+ 1 - x_size) && mod(ID_number_array(j), x_size+1) ~= 0 && ID_number_array(j) ~= 1
                p2 = coordinates(ID_number_array(j) + x_size - 1,:);
                C(ID_number_array(j), ID_number_array(j) + x_size - 1) = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));
                C(ID_number_array(j) + x_size - 1, ID_number_array(j)) = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));

            end
        end
        
    end
end

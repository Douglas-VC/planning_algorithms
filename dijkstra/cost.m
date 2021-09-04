function c = cost(p1,p2)
    global function_h;
    global function_c;
    global metric;    

    if strcmp(metric,'distance')
        c = abs(sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2));
        
    elseif strcmp(metric,'height')
        d = sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2);
        f = (function_c(p1) + function_c((p1+p2)/2) + function_c(p2))/3;
        c = abs(f*d);
        
    elseif strcmp(metric,'slope')
        c = abs((function_c(p1)+function_c(p2)+function_c((p1+p2)/2))/3 + max([function_c(p1),function_c(p2),function_c((p1+p2)/2)])/3);
        
    elseif strcmp(metric,'combined')
        d = sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 + (function_h(p1) - function_h(p2))^2);
        f = (function_c(p1) + function_c((p1+p2)/2) + function_c(p2))/3;
        c = abs(f*d + d);
    end    

%     costs = [];
%     costs(1,1) = d;
%     costs(2,1) = f*d;
%     costs(3,1) = total_angle;
%     costs(4,1) = f*d + d;

%     global all_costs;     
%     all_costs = [all_costs costs];

end
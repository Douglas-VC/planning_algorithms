
function f = f_cost(p)
    global function_h    
    
    f_0 = function_h(p(1),p(2));
    f_Mx = function_h(p(1)+1e-4,p(2));
    f_My = function_h(p(1),p(2)+1e-4);
    grad = ([f_My f_Mx]-f_0)/(1e-4);
    f = norm(grad);

end
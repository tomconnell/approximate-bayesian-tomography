function [logsmoothness] = smoothness(parameters,l,w)

total_sum = 0;

for ii = 1:length(parameters)
    
    holder_sum = 0;
    norm = 0;
    
    for jj = 1:4
        
        if jj == 1
            
            move = ii + 1;
            
            if move >= 1 && move <= l*w && ceil(move/l) == ceil(ii/l)
                norm = norm + 1;
                holder_sum = holder_sum +...
                    (parameters(move)-parameters(ii))^2;
            end
            
        end
        
        if jj == 2
            
            move = ii - 1;
            
            if move >= 1 && move <= l*w && ceil(move/l) == ceil(ii/l)
                norm = norm + 1;
                holder_sum = holder_sum +...
                    (parameters(move)-parameters(ii))^2;
            end
            
        end
        
        if jj == 3
            
            move = ii + l;
            
            if move >= 1 && move <= l*w
                
                norm = norm + 1;
                holder_sum = holder_sum +...
                    (parameters(move)-parameters(ii))^2;
            end
            
        end
        
        if jj == 4
            
            move = ii - l;
            
            if move >= 1 && move <= l*w
                
                norm = norm + 1;
                holder_sum = holder_sum +...
                    (parameters(move)-parameters(ii))^2;
            end
             
        end
        
    end
    
    total_sum = holder_sum/norm;
    
end

logsmoothness = total_sum;
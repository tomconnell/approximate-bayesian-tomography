function [logsmoothness] = smoothness(parameters)

total_sum = 0;

for ii = 1:length(parameters)
    
    holder_sum = 0;
    norm = 0;
    
    for jj = 1:4
        
        if jj == 1
            
            move = ii + 1;
            
            if move >= 1 && move <= 32 && ceil(move/4) == ceil(ii/4)
                norm = norm + 1;
                holder_sum = holder_sum +...
                    (parameters(move)-parameters(ii))^2;
            end
            
        end
        
        if jj == 2
            
            move = ii - 1;
            
            if move >= 1 && move <= 32 && ceil(move/4) == ceil(ii/4)
                norm = norm + 1;
                holder_sum = holder_sum +...
                    (parameters(move)-parameters(ii))^2;
            end
            
        end
        
        if jj == 3
            
            move = ii + 4;
            
            if move >= 1 && move <= 32
                
                norm = norm + 1;
                holder_sum = holder_sum +...
                    (parameters(move)-parameters(ii))^2;
            end
            
        end
        
        if jj == 4
            
            move = ii - 4;
            
            if move >= 1 && move <= 32
                
                norm = norm + 1;
                holder_sum = holder_sum +...
                    (parameters(move)-parameters(ii))^2;
            end
             
        end
        
    end
    
    total_sum = holder_sum/norm;
    
end

logsmoothness = total_sum;
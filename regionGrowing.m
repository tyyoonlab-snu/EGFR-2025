function [region, finalImage, raw_data] = regionGrowing(image, x_seed, y_seed, T, raw_image)
    [rows, cols] = size(image);
    visited = zeros(rows, cols);
    region = zeros(rows, cols);
    stack = [x_seed, y_seed];
    seed_value = image(x_seed, y_seed);
    
    while ~isempty(stack)
        % Pop a point from the stack
        point = stack(1, :);
        stack(1, :) = [];
        
        x = point(1);
        y = point(2);
        
        % Check if the point has been visited
        if visited(x, y)
            continue;
        end
        
        visited(x, y) = 1;
        
        % Check if the point is similar to the seed
        dif =abs(double(image(x, y)) - double(seed_value));
        if dif <= T && double(image(x,y))~=0
            region(x, y) = 1;
            
            % Add neighbors to stack
            for dx = -1:1
                for dy = -1:1
                    x_new = x + dx;
                    y_new = y + dy;
                    
                    if x_new >= 1 && x_new <= rows && y_new >= 1 && y_new <= cols
                        stack = [stack; x_new, y_new];
                    end
                end
            end
        end
    end
    
    finalImage = double(region) .* double(image);
    raw_data = double(region) .* double(raw_image);
    
end
function [ disparity ] = stereo_computation( Il, Ir, min_dis, max_dis, wsize, cost_function )
% Input parameters are 5:
%     - Il: left image
%     - Ir: right image
%     - min_dis: minimum disparity
%     - max_dis: maximum disparity
%     - wsize: window size (e.g. a value of 3 indicates a 3x3 window)
%     - cost_function: matching cost (the user may able to choose between SSD and NCC costs)

[h,w] = size(Il);
disparity=zeros(h-wsize+1,w-wsize+1);
i_dis=1;
for i=1:h-wsize+1
    j_dis = 1;
    for j=1:h-wsize+1
        % Get current window
        reference_window = Il(i:i+wsize-1,j:j+wsize-1);
        
        % Compute the matching cost
        if upper(cost_function)=='SSD'
            costs = zeros(w-wsize+1,1);
            for k=1:w-wsize+1
                comp_window = Ir(i:i+wsize-1,k:k+wsize-1);
                costs(k) = sum(sum((reference_window-comp_window).^2));
            end
        else
            error(strcat(cost_function,' not a valid matching cost name.'));
        end
        
        % Find minimum
        [~,match_pos] = min(costs);
        match_pos = match_pos + floor(wsize/2);

        % Compute disparity
        current_disparity = i-match_pos;
        if current_disparity<min_dis, current_disparity=min_dis; end
        if current_disparity>max_dis, current_disparity=max_dis; end
        
        disparity(i_dis,j_dis) = current_disparity;
        
        j_dis=j_dis+1;
    end
    i_dis = i_dis+1;
end

end


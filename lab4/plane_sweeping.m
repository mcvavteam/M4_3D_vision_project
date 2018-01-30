function [ I_depth ] = plane_sweeping( I,P1,P2,window_size,threshold,matching_function,plot_im_matching )

weights = ones(window_size)/window_size.^2;
[rows,cols] = size(I{1});
corners = [1, cols, 1, rows];
I_depth = inf(rows,cols);
best_matching = inf(rows,cols);
pad = floor(window_size/2);
I_1 = padarray(I{1},[pad,pad],'replicate');
for depth=1:15
    % Compute the fronto-parallel plane at the corresponding depth
    fronto_parallel_plane = P2(3,:)-[0,0,0,depth];
    
    % Compute the image to image homography given a plane
    A = inv([P2; fronto_parallel_plane]);
    H = P1*A(:,1:3);
    
    % Warp the second image according to the homography
    I_2projected = apply_H_v2(I{2}, H, corners);
    if plot_im_matching
        auxI = repmat(I{1},[1,1,2]);
        auxI(:,:,3) = I_2projected;
        %auxI(auxI==Inf)=0;
        imshow(auxI);
        title(['d_', num2str(depth)]);
        pause
    end
    
    % Compute the matching score (SSD, NCC, ...)
    I_2projected = padarray(I_2projected,[pad,pad],'replicate');
    I_2projected(I_2projected==0)=Inf;
    matching = zeros(rows,cols);
    for i=1+pad:rows+pad
        cur_I1row = zeros(window_size,window_size,cols);
        cur_I2row = zeros(window_size,window_size,cols);
        for wj=1:window_size
            cur_I1row(:,wj,:) = permute(I_1(i-pad:i+pad,wj:cols+(2*pad)-window_size+wj),[1,3,2]);
            cur_I2row(:,wj,:) = permute(I_2projected(i-pad:i+pad,wj:cols+(2*pad)-window_size+wj),[1,3,2]);
        end
        
        matching(i-pad,:) = permute(sum(sum(bsxfun(@times,weights,...
            (cur_I1row-cur_I2row).^2),1),2),[1,3,2]);
    end
    
    % For each pixel keep the depth with best score
    I_depth(best_matching > matching & matching < threshold) = depth;
    best_matching(best_matching > matching & matching < threshold) = matching(best_matching > matching & matching < threshold);
end
I_depth(I_depth==Inf)=0;
end


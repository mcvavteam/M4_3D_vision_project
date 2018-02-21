function [ disparity,costs ] = stereo_computation_bp( Il, Ir, min_dis, max_dis, wsize, cost_function, bilateral)
% Input parameters are 5:
%     - Il: left image
%     - Ir: right image
%     - min_dis: minimum disparity
%     - max_dis: maximum disparity
%     - wsize: window size (e.g. a value of 3 indicates a 3x3 window)
%     - cost_function: matching cost (the user may able to choose between SSD and NCC costs)

    [h,w] = size(Il);
    disparity=zeros(h-wsize,w-wsize);
    costs=Inf((h-wsize)*(w-wsize),33);
    
    w_step = floor(wsize/2);
    count = 1;
    Il = padarray(Il, [w_step w_step]);
    Ir = padarray(Ir, [w_step w_step]);

    for i=1+w_step:h-w_step
        for j=1+w_step:w-wsize
            windowl = double(Il(i-w_step:i+w_step, j-w_step:j+w_step));

            kmin = max(1+w_step, j-max_dis);
            kmax = min(w-w_step, j+max_dis);
                        
            [wh,ww] = size(windowl);
            weight = repmat(1/(wh*ww),wh,ww);
%                 min_error = Inf;
%                 for k = kmin:kmax
%                     windowr = double(Ir(i-w_step:i+w_step, k-w_step:k+w_step));
%                     err_k = sum(sum(weight.*(windowl-windowr).^2));
%                     if err_k < min_error
%                         min_error = err_k;
%                         kbest = k;
%                     end
%                 end

            err = zeros(wsize.^2,kmax-kmin+1);
            for k = kmin:kmax
                windowr = double(Ir(i-w_step:i+w_step, k-w_step:k+w_step));
                windowrp(:,k-kmin+1) = windowr(:);
            end
            err = bsxfun(@minus,windowrp,windowl(:));
            err = bsxfun(@times,err.^2,weight(:));
            err = sum(err,1);
%             [~, kbest] = min(err);
%             kbest = kbest+kmin-1;
%             if(j-max_dis < 1+w_step)
%                 err = [Inf(1,1+w_step-(j-max_dis)),err];
%             end
%             if(w-w_step < j+max_dis)
%                 err = [err, Inf(1,j+max_dis-w-w_step)];
%             end
            if(j-max_dis < 1+w_step)
                err = [Inf(1,33-size(err,2)),err];
            end
            if(w-w_step < j+max_dis)
                err = [err, Inf(1,33-size(err,2))];
            end
            [ ~, index] = min(err);
            costs(((w-wsize)*(i-1))+j,:) = err; 
            disparity(i,j) = index;
%             disparity(i-w_step, j-w_step) = abs(j-kbest);
        end
        disp(i)
    end
end


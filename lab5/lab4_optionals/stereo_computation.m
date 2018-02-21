function [ disparity ] = stereo_computation( Il, Ir, min_dis, max_dis, wsize, cost_function, bilateral)
% Input parameters are 5:
%     - Il: left image
%     - Ir: right image
%     - min_dis: minimum disparity
%     - max_dis: maximum disparity
%     - wsize: window size (e.g. a value of 3 indicates a 3x3 window)
%     - cost_function: matching cost (the user may able to choose between SSD and NCC costs)

    [h,w] = size(Il);
    disparity=zeros(h-wsize,w-wsize);
    
    w_step = floor(wsize/2);
    
    for i=1+w_step:h-w_step
        for j=1+w_step:w-wsize
            windowl = double(Il(i-w_step:i+w_step, j-w_step:j+w_step));

            kmin = max(1+w_step, j-max_dis);
            kmax = min(w-w_step, j+max_dis);
                        
            [wh,ww] = size(windowl);
            weight = repmat(1/(wh*ww),wh,ww);
                 
            % Bilateral Weights
            if(bilateral)
                gammac = 5;
                gammap = w_step;
                pixel_center = double(Il(i,j));
                Ac = windowl-pixel_center;
                [grid1, grid2] = meshgrid(i-w_step:i+w_step, j-w_step:j+w_step);
                Ag = sqrt(((grid1-i).^2 +(grid2-j).^2));

                weight = exp(-Ac/gammac).*exp(-Ag/gammap);
            end
            
            if strcmpi(cost_function,'SSD')                      
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
                [~, kbest] = min(err);
                kbest = kbest+kmin-1;

                
            elseif strcmpi(cost_function,'NCC')  
                min_error = -Inf;
                for k = kmin:kmax
                    windowr = double(Ir(i-w_step:i+w_step, k-w_step:k+w_step));
                    
                    meanl = sum(sum(windowl.*weight));
                    meanr = sum(sum(windowr.*weight));                   
                    sigl = sqrt(sum(sum(weight.*((windowl-meanl).^2))));
                    sigr = sqrt(sum(sum(weight.*((windowr-meanr).^2))));
                    
                    err_k = sum(sum(weight.*(windowl-meanl).*(windowr-meanr)))/(sigl*sigr);

                    if err_k > min_error
                        min_error = err_k;
                        kbest = k;
                    end
                end
            else
                error(strcat(cost_function,' not a valid matching cost name.'));
            end
            disparity(i-w_step, j-w_step) = (j-kbest);
        end
    end
end


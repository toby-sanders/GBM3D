function tau = getTau(levels,mode)

% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020

% empirically tuned threshold values for wavelet shrinkage 

% traditionally tau = 3.3 - .3*(i-1) worked reasonably well
tau = zeros(levels,7);
switch mode
    case 1
        for i = 1:levels
            tau(i,:) = 3.3- 0.3*(i-1);
        end
    case 2        
        for j = 1:7
            for i = 1:levels
                if j<4
                    tau(i,j) = 3.3 - .3*(i-1);
                elseif j<7
                    tau(i,j) = 3.3 - .3*(i-1) + .5;
                else
                    tau(i,j) = 3.3 - .3*(i-1) + 1;
                end
            end
        end
    case 3
        for j = 1:7
            for i = 1:levels
                if j<4
                    tau(i,j) = 3.3 - .3*(i-1);
                elseif j<7
                    tau(i,j) = 3.3 - .3*(i-1) + .25;
                else
                    tau(i,j) = 3.3 - .3*(i-1) + 1.5;
                end
            end
        end
    case 4
        for j = 1:7
            for i = 1:levels
                if j<4
                    tau(i,j) = 2.9 ;
                elseif j<7
                    tau(i,j) = 2.9  + .5;
                else
                    tau(i,j) = 2.9  + 1;
                end
            end
        end        
end


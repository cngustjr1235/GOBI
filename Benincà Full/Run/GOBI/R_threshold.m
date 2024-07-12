function [R_tmp] = R_threshold(R, thres)
R_tmp = R;
for i = 1:length(R)
    if R(i) < thres
        R_tmp(i) = 0;
    else
        R_tmp(i) = 1;
        %L_tmp(i,j) = L_tmp(i,j);
    end
end
end


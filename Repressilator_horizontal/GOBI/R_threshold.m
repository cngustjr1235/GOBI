function [R_tmp] = R_threshold(R, thres)
R_tmp = R;
for i = 1:length(R)
    if isnan(R(i))
        R_tmp(i) = 0;
    elseif R(i) < thres
        R_tmp(i) = 0;
    else
        R_tmp(i) = 1;
        %R_tmp(i,j) = R_tmp(i,j);
    end
end
end


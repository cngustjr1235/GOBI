function [score_list, t_1, t_2] = RDS_ps_dim3(X, Y, Z, W, t, time_interval)

f = gradient(W, time_interval);
score_list = zeros(length(t),length(t),8);
t_1 = zeros(length(t),length(t));
t_2 = zeros(length(t),length(t));

for t_ori = 1:length(t)
    for t_prime = 1:length(t)
        t_1(t_ori, t_prime) = t(t_ori);
        t_2(t_ori, t_prime) = t(t_prime);
        
        X_d = X(t_ori) - X(t_prime);
        Y_d = Y(t_ori) - Y(t_prime);
        Z_d = Z(t_ori) - Z(t_prime);
        W_d = W(t_ori) - W(t_prime);
        f_d = f(t_ori) - f(t_prime);
        score_tmp = X_d * Y_d * Z_d * W_d * f_d;
        
        if     X_d >= 0 && Y_d >= 0 && Z_d >= 0 && W_d >= 0
            score_list(t_ori, t_prime,1) = score_tmp;
        elseif X_d >= 0 && Y_d >= 0 && Z_d < 0 && W_d >= 0
            score_list(t_ori, t_prime,2) = -score_tmp;
        elseif X_d >= 0 && Y_d < 0 && Z_d >= 0 && W_d >= 0
            score_list(t_ori, t_prime,3) = -score_tmp;
        elseif X_d >= 0 && Y_d < 0 && Z_d < 0 && W_d >= 0
            score_list(t_ori, t_prime,4) = score_tmp;
        elseif X_d < 0 && Y_d >= 0 && Z_d >= 0 && W_d >= 0
            score_list(t_ori, t_prime,5) = -score_tmp;
        elseif X_d < 0 && Y_d >= 0 && Z_d < 0 && W_d >= 0
            score_list(t_ori, t_prime,6) = score_tmp;
        elseif X_d < 0 && Y_d < 0 && Z_d >= 0 && W_d >= 0
            score_list(t_ori, t_prime,7) = score_tmp;
        elseif X_d < 0 && Y_d < 0 && Z_d < 0 && W_d >= 0
            score_list(t_ori, t_prime,8) = -score_tmp;
        end
    end
end
end
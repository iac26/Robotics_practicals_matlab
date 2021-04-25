function [x_hat,P_hat] = BARO_update_const_a(p,x_tilde,P_tilde,R_baro,g)

    h = x_tilde(1);
    p0 = x_tilde(4);
    k = x_tilde(5);
    h0 = x_tilde(6);
    
    p_est = p0*exp(k*g*(h0-h));

    H = [-k*g*p_est, 0, 0, p_est/p0, g*(h0-h)*p_est,  k*g*p_est];
    
    K = P_tilde*H'*1/(H*P_tilde*H' + R_baro);
    
    x_hat = x_tilde + K*(p - p_est);
    
    P_hat = (eye(6) - K*H)*P_tilde;
    
end




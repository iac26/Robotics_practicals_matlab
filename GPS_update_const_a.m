function [x_hat,P_hat] = GPS_update_const_a(z,x_tilde,P_tilde,R_gps)
    
    H = [1, 0, 0, 0, 0, 0];
    
    K = P_tilde*H'/(H*P_tilde*H' + R_gps);
    
    x_hat = x_tilde + K*(z - H*x_tilde);
    
    P_hat = (eye(6) - K*H)*P_tilde;


end


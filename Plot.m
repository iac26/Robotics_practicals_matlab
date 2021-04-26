figure('Renderer', 'painters', 'Position', [10 10 900 600]);
subplot(221);
hold on

top1 = X(1,:)+3*sqrt(P_diag(1,:));
bot1 = X(1,:)-3*sqrt(P_diag(1,:));

plot(X(1,:), 'b');
plot(top1, 'r-');
plot(bot1, 'r-');
title("h");

subplot(222);
hold on
top4 = X(4,:)+3*sqrt(P_diag(4,:));
bot4 = X(4,:)-3*sqrt(P_diag(4,:));

plot(X(4,:), 'b');
plot(top4, 'r-');
plot(bot4, 'r-');
title("p_0");

subplot(223);
hold on
top5 = X(5,:)+3*sqrt(P_diag(5,:));
bot5 = X(5,:)-3*sqrt(P_diag(5,:));

plot(X(5,:), 'b');
plot(top5, 'r-');
plot(bot5, 'r-');
title("K");

subplot(224);
hold on
top6 = X(6,:)+3*sqrt(P_diag(6,:));
bot6 = X(6,:)-3*sqrt(P_diag(6,:));

plot(X(6,:), 'b');
plot(top6, 'r-');
plot(bot6, 'r-');
title("h_0");

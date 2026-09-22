// MUX recebendo as saidas dos ganhos 1/R1 e 1/R2
plot(A.time, A.values(:,1), 'b', A.time, A.values(:,2), 'r');
xlabel('t (s)'); ylabel('vazao (m^3/s)');
legend('q1','q2'); xgrid();
disp(A.values($,:));   // 0.1   0.1

// MUX recebendo (q - q2) e (C1*dh1 + C2*dh2)
plot(A.time, A.values(:,1), 'b', A.time, A.values(:,2), 'r--');
xlabel('t (s)'); ylabel('vazao (m^3/s)');
legend('q - q2', 'C1*dh1 + C2*dh2'); xgrid();
disp(max(abs(A.values(:,1) - A.values(:,2))));   // 8.743e-17

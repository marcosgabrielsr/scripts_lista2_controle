// MUX recebendo as saidas dos dois integradores
plot(A.time, A.values(:,1), 'k', A.time, A.values(:,2), 'b');
xlabel('t (s)'); ylabel('nivel (m)');
legend('h1','h2'); xgrid();
disp(A.values($,:));   // 0.5   0.3

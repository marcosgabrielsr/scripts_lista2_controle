// =====================================================================
// Exercicio 3 - Simulacao da planta de 6a ordem e reducoes de ordem
// Os trechos abaixo sao executados em sequencia, compartilhando as
// variaveis definidas no primeiro deles.
// =====================================================================

// ---------------------------------------------------------------------
// Definicao da planta e dos modelos reduzidos
// ---------------------------------------------------------------------
// Planta completa - 6a ordem
s  = %s;
n  = s + 6;
d6 = s^6 + 11*s^5 + 52*s^4 + 142*s^3 + 241*s^2 + 231*s + 90;
g6 = syslin('c', n, d6);
disp(roots(d6));   // -1, -1+-2i, -2, -3, -3
disp(roots(n));    // -6

// Modelos reduzidos - ganho estatico preservado
d5 =  3*(s^5 + 8*s^4 + 28*s^3 + 58*s^2 + 67*s + 30);
d4 =  9*(s^4 + 5*s^3 + 13*s^2 + 19*s + 10);
d3 = 18*(s^3 + 3*s^2 + 7*s + 5);
g5 = syslin('c', n, d5);
g4 = syslin('c', n, d4);
g3 = syslin('c', n, d3);

// Mapa de polos e zeros
p = roots(d6);  z = roots(n);
scf(0); plot(real(p), imag(p), 'bx'); plot(real(z), imag(z), 'ro'); xgrid();

// ---------------------------------------------------------------------
// Simulacao do modelo completo - item (a)
// ---------------------------------------------------------------------
// Entradas u = t^i, i = 0,1,2
t = 0 : .005 : 12;  N = length(t);
U = [ones(1,N); t; t.^2];

// Assintotas previstas pelo teorema do valor final
K0 = 6/90;  K1 = -0.16;  K2 = 0.464296;
A  = [K0*ones(1,N); K0*t + K1; K0*t.^2 + 2*K1*t + K2];

for i = 1:3
    y6 = csim(U(i,:), t, g6);
    scf(i); plot(t, y6, 'b', t, A(i,:), 'k--'); xgrid();
    legend('modelo completo', 'assintota');
end

// ---------------------------------------------------------------------
// Comparacao com os modelos reduzidos - item (c)
// ---------------------------------------------------------------------
for i = 1:3
    y6 = csim(U(i,:), t, g6);
    y5 = csim(U(i,:), t, g5);
    y4 = csim(U(i,:), t, g4);
    y3 = csim(U(i,:), t, g3);

    scf(3+i);
    subplot(2,1,1);
    plot(t, y6, 'k', t, y5, 'b', t, y4, 'g', t, y3, 'r'); xgrid();
    legend('completa', '5a ordem', '4a ordem', '3a ordem');
    subplot(2,1,2);
    plot(t, y6-y5, 'b', t, y6-y4, 'g', t, y6-y3, 'r'); xgrid();
    legend('5a ordem', '4a ordem', '3a ordem');

    // max|e(t)| e e(t_final) de cada reducao
    disp([max(abs(y6-y5)), y6($)-y5($)]);
    disp([max(abs(y6-y4)), y6($)-y4($)]);
    disp([max(abs(y6-y3)), y6($)-y3($)]);
end

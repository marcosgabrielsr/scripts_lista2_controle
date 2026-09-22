// =====================================================================
// Exercicio 4.c) - Reducao de ordem por g(0) e g'(0)
// =====================================================================
clear; clc;
s = %s;
mkdir('Imagens');

deng = poly([16875 27000 15075 5720 817 48 1], 's', 'coeff');
g    = syslin('c', 1, deng);

g0  = horner(g, 0);                       // 5.9259e-5
dg0 = horner(derivat(g), 0);              // -9.4815e-5

// gr(s) = (b1 s + b0)/(s^3+3s^2+7s+5):  gr(0) = b0/5 ; gr'(0) = (5b1-7b0)/25
b0 = 5*g0;                                // 2.9630e-4
b1 = (25*dg0 + 7*b0)/5;                   // -5.9259e-5
gr = syslin('c', b1*s + b0, (s + 1)*(s^2 + 2*s + 5));
disp(b0, b1);
disp(roots(b1*s + b0));                   // zero em +5 (fase nao minima)
disp(horner(gr, 0) - g0, horner(derivat(gr), 0) - dg0);   // ~0

t  = 0:0.01:8;
yg = csim('step', t, g);
yr = csim('step', t, gr);
disp(max(abs(yg - yr)));                  // ~1e-6

scf(6); clf();
plot(t, yg, 'b'); plot(t, yr, 'r--'); xgrid();
legend(['planta completa'; 'planta reduzida'], 4);
xtitle('Item (c) - reducao de ordem', 't [s]', 'y(t)');
xs2png(6, 'Imagens/fig_ex4c_reducao.png');

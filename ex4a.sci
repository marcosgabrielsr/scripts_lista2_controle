// =====================================================================
// Exercicio 4.a)
// Projeto de compensador via equacao diofantina para
//     g(s) = 1 / (s^6 + 48s^5 + 817s^4 + 5720s^3 + 15075s^2 + 27000s + 16875)
// com: erro estatico de posicao nulo, dominancia de 2a ordem,
//      beta < 40 graus e sigma > 3.
// =====================================================================
clear; clc;

s = %s;

// -------------------------------------------------------- 1. Planta
coef_g = [16875, 27000, 15075, 5720, 817, 48, 1];
deng   = poly(coef_g, 's', 'coeff');
g      = syslin('c', 1, deng);
disp(roots(deng));      // -1, -1+-2i e -15 (triplo)

// ------------------------------------------------ 2. Planta aumentada
// planta tipo 0 (deng(0) = 16875): erro de posicao so zera com
// integrador no compensador -> g~ = g/s ; d_g~ = s*deng, grau ng = 7
dengt = s*deng;
disp(degree(dengt)); // 7

// ---------------------------------------- 3. Polos desejados de MF
sig = 3.5;  bet = 35*%pi/180;
wd  = sig*tan(bet);
p1  = -sig + wd*%i;  p2 = -sig - wd*%i;

// 10 polos rapidos distintos (>=5x mais afastados que sigma)
prapido = -(18:27);
polos = [p1, p2, prapido];
q = poly(polos, 's', 'roots');

disp(wd); disp(cos(bet));   // wd = 2.4507 ; zeta = 0.8192
disp(degree(q));            // 12

// ------------------ 4. Dimensoes da matriz de Sylvester (S quadrada)
// mc = ng - 1 = 6 ; propriedade de c exige nc = 5 ; L = nc+ng+1 = 13
nc = 5;  mc = 6;  L = nc + 7 + 1;

// ------------- 5. Solucao por divisao polinomial (evita mau condic.)
[numc, denc] = pdiv(q, dengt);   // q = dengt*denc + numc
disp(denc);  // d_c~ (grau 5)
disp(numc);  // n_c~ (grau 6)

// ------------------------------- 6. Compensador e malha fechada
c    = syslin('c', numc, s*denc);
FTMF = (c*g) /. 1;

disp(roots(denc));      // polos do compensador (estavel)
disp(roots(numc));      // zeros de malha fechada
disp(horner(FTMF, 0));  // ganho DC = 1 -> erro de posicao nulo

// ------------------------------------------- 7. Resposta ao degrau
t     = 0:0.001:2.5;
y_dio = csim('step', t, FTMF);

wn = sig/cos(bet);
s2 = syslin('c', wn^2, s^2 + 2*cos(bet)*wn*s + wn^2);
y_2a = csim('step', t, s2);

disp(max(y_dio) - 1);   // sobressinal real     (~0.895)
disp(max(y_2a)  - 1);   // sobressinal previsto (~0.0113)

mkdir('Imagens');
scf(11); clf(); plot(t, y_dio, 'b'); plot(t, y_2a, 'k--'); xgrid();
legend(['diofantina (malha fechada completa)'; '2a ordem dos polos dominantes']);
xtitle('4a - resposta ao degrau', 't [s]', 'y(t)');
xs2png(11, 'Imagens/fig_ex4a_resposta.png');

// -------------------- Zona de desempenho garantido (Omega) no plano s
b40 = 40*%pi/180;  sgm = 3;  Lp = 30;
scf(12); clf();
a = gca(); a.data_bounds = [-Lp, -12; 2, 12]; a.tight_limits = "on"; a.box = "on";
xv = [-sgm, -Lp, -Lp, -sgm];
yv = [ sgm*tan(b40), Lp*tan(b40), -Lp*tan(b40), -sgm*tan(b40)];
xfpoly(xv, yv);
e = gce(); e.background = color(220,235,250); e.clip_state = "clipgrf";
plot(real(polos), imag(polos), 'bx');
z = roots(numc);
plot(real(z), imag(z), 'ro');
plot([-sig -sig], [wd -wd], 'k*');
a.data_bounds = [-Lp, -12; 2, 12];
xgrid(); xtitle('Plano s - zona de desempenho garantido', 'Re', 'Im');
xs2png(12, 'Imagens/fig_ex4a_omega.png');

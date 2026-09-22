// =====================================================================
// Exercicio 4.e) - Comparacao dos projetos na planta completa
// =====================================================================
clear; clc;
s = %s;
mkdir('Imagens');

function a = angd(x)
    a = atan(imag(x), real(x))*180/%pi;
endfunction
function a = wrap(x)
    a = pmodulo(x + 180, 360) - 180;
endfunction
function ts = tacom(t, y)                 // tempo de acomodacao (2%)
    idx = find(abs(y - 1) > 0.02);
    ts  = t(max(idx));
endfunction

deng = poly([16875 27000 15075 5720 817 48 1], 's', 'coeff');
g    = syslin('c', 1, deng);
nr   = (5 - s)/16875;
sig  = 3.5;  bet = 35;
sd   = -sig + %i*sig*tand(bet);

// ---- (a) diofantina
q = poly([sd, conj(sd), -(18:27)], 's', 'roots');
[numc_a, denc_a] = pdiv(q, s*deng);
c_a = syslin('c', numc_a, s*denc_a);

// ---- (b) anulacao + avanco, projetado na planta completa
gc  = syslin('c', 1, s*(s + 81)^2*(s + 15)^3);
phi = wrap(180 - angd(horner(gc, sd)));
zb  = -14;  pb = real(sd) - imag(sd)/tand(angd(sd - zb) - phi);
kb  = 1/abs(horner(gc, sd)*(sd - zb)/(sd - pb));
c_b = syslin('c', kb*(s + 1)*(s^2 + 2*s + 5)*(s - zb), s*(s + 81)^2*(s - pb));

// ---- (d) anulacao + atraso, projetado na planta reduzida
gcr = syslin('c', nr, s*(s + 81)^2);
phi = wrap(180 - angd(horner(gcr, sd)));
zd  = -14;  pd = real(sd) - imag(sd)/tand(angd(sd - zd) - phi);
kd  = 1/abs(horner(gcr, sd)*(sd - zd)/(sd - pd));
c_d = syslin('c', kd*(s + 1)*(s^2 + 2*s + 5)*(s - zd), s*(s + 81)^2*(s - pd));

// ---- malhas fechadas na planta COMPLETA
Fa = (c_a*g) /. 1;   Fb = (c_b*g) /. 1;   Fd = (c_d*g) /. 1;

t  = 0:0.001:4;
ya = csim('step', t, Fa);  yb = csim('step', t, Fb);  yd = csim('step', t, Fd);
mprintf("(a) Mp=%.2f%%  ts=%.2f s\n", (max(ya)-1)*100, tacom(t, ya));
mprintf("(b) Mp=%.2f%%  ts=%.2f s\n", (max(yb)-1)*100, tacom(t, yb));
mprintf("(d) Mp=%.2f%%  ts=%.2f s\n", (max(yd)-1)*100, tacom(t, yd));
disp(roots(s*(s + 81)^2*(s - pd)*(s + 15)^3 + kd*(s - zd)));  // par -2.26+-2.42i

scf(10); clf();
plot(t, ya, 'k'); plot(t, yb, 'b'); plot(t, yd, 'r--'); xgrid();
set(gca().children.children, "thickness", 2);
legend(['(a) diofantina'; '(b) LR - planta completa'; ..
        '(d) LR - planta reduzida'], 1);
xtitle('Item (e) - comparacao na planta completa', 't [s]', 'y(t)');
xs2png(10, 'Imagens/fig_ex4e_comparacao.png');

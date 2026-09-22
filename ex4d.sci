// =====================================================================
// Exercicio 4.d) - Lugar das raizes na planta reduzida
//   gr(s) = (5 - s) / (16875 (s+1)(s^2+2s+5))
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

deng = poly([16875 27000 15075 5720 817 48 1], 's', 'coeff');
g    = syslin('c', 1, deng);
nr   = (5 - s)/16875;                       // numerador de gr
dr   = (s + 1)*(s^2 + 2*s + 5);
gr   = syslin('c', nr, dr);
gri  = syslin('c', nr, s*dr);               // gr(s)/s

sig = 3.5;  bet = 35;
sd  = -sig + %i*sig*tand(bet);

// ------------------------------------------------ Parte 1 (LR de 0 graus)
nn = s - 5;  dd = s*dr;
bq = derivat(nn)*dd - nn*derivat(dd);       // -3s^4+14s^3+38s^2+70s+25
cand = roots(bq);
for k = 1:size(cand, '*')
    mprintf("candidato %s   K = %s\n", string(cand(k)), ..
            string(-horner(dd, cand(k))/horner(nr, cand(k))));
end
disp(routh_t(gri, poly(0, 'k')));
Kp = kpure(gri);  disp(min(Kp(Kp > 0)));    // 23521.7

scf(7); clf();
evans(gri, 40000);                          // gr ja contem o sinal negativo
xtitle('LR de gr(s)/s', 'Re(s)', 'Im(s)');
xs2png(7, 'Imagens/fig_ex4d_lr_p.png');

// ------------------------------------------------ Parte 2
phi2 = wrap(180 - angd(horner(gri, sd)));   // 45.76
z2   = -3.5;
p2   = real(sd) - imag(sd)/tand(angd(sd - z2) - phi2);   // -6.0163
kc2  = 1/abs(horner(gri, sd)*(sd - z2)/(sd - p2));      // 530341
h2   = syslin('c', (s - z2)*nr, s*(s - p2)*dr);
Kp    = kpure(h2);  Klim2 = min(Kp(Kp > 0)); // 43771.5
mprintf("phi=%.2f  p=%.4f  kc=%.1f  Klim=%.1f  kc/Klim=%.2f\n", ..
        phi2, p2, kc2, Klim2, kc2/Klim2);
disp(roots(s*(s - p2)*dr + kc2*(s - z2)*nr));   // par +1.41 +- 2.07i

scf(8); clf();
evans(h2, 7e5);
xtitle('LR da planta reduzida com integrador + avanco', 'Re(s)', 'Im(s)');
xs2png(8, 'Imagens/fig_ex4d_lr_avanco.png');

// ------------------------------------------------ Parte 3
gcr  = syslin('c', nr, s*(s + 81)^2);       // apos o cancelamento
phi3 = wrap(180 - angd(horner(gcr, sd)));   // -15.29
z3   = -4*sig;                              // -14
p3   = real(sd) - imag(sd)/tand(angd(sd - z3) - phi3);   // -8.0265 (atraso)
kc3  = 1/abs(horner(gcr, sd)*(sd - z3)/(sd - p3));      // 2.3394e7
mprintf("phi=%.2f  z=%.2f  p=%.4f  kc=%.4e\n", phi3, z3, p3, kc3);

c_d     = syslin('c', kc3*(s + 1)*(s^2 + 2*s + 5)*(s - z3), ..
                       s*(s + 81)^2*(s - p3));
FTMF_dr = (c_d*gr) /. 1;                    // no modelo reduzido
FTMF_df = (c_d*g)  /. 1;                    // na planta completa
disp(roots(s*(s + 81)^2*(s - p3) + kc3*(s - z3)*nr));   // reduzida
disp(roots(s*(s + 81)^2*(s - p3)*(s + 15)^3 + kc3*(s - z3)));  // completa

t    = 0:0.001:4;
y_dr = csim('step', t, FTMF_dr);
y_df = csim('step', t, FTMF_df);
mprintf("reduzida: Mp=%.2f%%  subsinal=%.1f%%\n", (max(y_dr)-1)*100, -min(y_dr)*100);
mprintf("completa: Mp=%.2f%%\n", (max(y_df)-1)*100);

scf(9); clf();
plot(t, y_dr, 'r'); plot(t, y_df, 'b--'); xgrid();
legend(['modelo reduzido'; 'planta completa'], 4);
xtitle('Item (d) - resposta ao degrau', 't [s]', 'y(t)');
xs2png(9, 'Imagens/fig_ex4d_resposta.png');

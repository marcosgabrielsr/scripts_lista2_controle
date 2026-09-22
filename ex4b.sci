// =====================================================================
// Exercicio 4.b) - Lugar das raizes na planta completa
//   Parte 1: LR de g(s)/s (proporcional + integrador)
//   Parte 2: integrador + avanco (falha, verificado por Routh)
//   Parte 3: anulacao dos polos com Re=-1 + avanco
// =====================================================================
clear; clc;
s = %s;
mkdir('Imagens');

function a = angd(x)          // angulo em graus
    a = atan(imag(x), real(x))*180/%pi;
endfunction
function a = wrap(x)          // leva o angulo para (-180, 180]
    a = pmodulo(x + 180, 360) - 180;
endfunction
function plot_omega(sg, bt, polos, zs, xmin, ymax)   // regiao Omega
    a  = gca();
    xv = [-sg, xmin, xmin, -sg];
    yv = [sg*tand(bt), -xmin*tand(bt), xmin*tand(bt), -sg*tand(bt)];
    a.data_bounds  = [xmin, -ymax; 2, ymax];
    a.tight_limits = "on";  a.box = "on";
    xfpoly(xv, yv); e = gce();
    e.background = color(220,235,250); e.foreground = color(120,120,120);
    e.clip_state = "clipgrf";
    plot(real(polos), imag(polos), 'bx');
    if ~isempty(zs) then plot(real(zs), imag(zs), 'ro'); end
    a.data_bounds  = [xmin, -ymax; 2, ymax];
    xgrid();
endfunction

deng = poly([16875 27000 15075 5720 817 48 1], 's', 'coeff');
g    = syslin('c', 1, deng);
gi   = syslin('c', 1, s*deng);                  // planta + integrador

sig = 3.5;  bet = 35;
sd  = -sig + %i*sig*tand(bet);                  // -3.5 + j2.4507

// ------------------------------------------------ Parte 1
disp(roots(s*deng));                            // 0, -1, -1+-2i, -15 (x3)
dd   = derivat(s*deng);                         // n'd - nd' = -d'(s)
cand = roots(dd);
for k = 1:size(cand, '*')
    mprintf("candidato %s   K = %s\n", string(cand(k)), ..
            string(-horner(s*deng, cand(k))));
end
disp(krac2(gi));                                // ganho no ponto de quebra (3280.78)
disp(routh_t(gi, poly(0, 'k')));                // arranjo de Routh em k
Kp = kpure(gi);  disp(min(Kp(Kp > 0)));         // cruzamento: 24204.74

scf(1); clf();
evans(gi, 60000);
a = gca(); a.data_bounds = [-18, -4; 2, 4]; a.tight_limits = "on";
xtitle('LR de g(s)/s', 'Re(s)', 'Im(s)');
xs2png(1, 'Imagens/fig_ex4b_lr_p.png');

// ------------------------------------------------ Parte 2
phi2 = wrap(180 - angd(horner(gi, sd)));        // 65.76
z2   = -3.5;
p2   = real(sd) - imag(sd)/tand(angd(sd - z2) - phi2);        // -8.9439
kc2  = 1/abs(horner(gi, sd)*(sd - z2)/(sd - p2));             // 768170
h2   = syslin('c', (s - z2), s*(s - p2)*deng);
Kp    = kpure(h2);  Klim2 = min(Kp(Kp > 0));    // 70520.04
mprintf("phi=%.2f  p=%.4f  kc=%.1f  Klim=%.2f  kc/Klim=%.2f\n", ..
        phi2, p2, kc2, Klim2, kc2/Klim2);
disp(roots(s*(s - p2)*deng + kc2*(s - z2)));    // par +1.28 +- 2.19i

scf(2); clf();
evans(h2, 1e6);
xtitle('LR com integrador + avanco (Parte 2)', 'Re(s)', 'Im(s)');
xs2png(2, 'Imagens/fig_ex4b_lr_avanco.png');

// ------------------------------------------------ Parte 3
gcanc = syslin('c', 1, s*(s + 81)^2*(s + 15)^3);   // apos o cancelamento
phi3  = wrap(180 - angd(horner(gcanc, sd)));        // 4.71
z3    = -4*sig;                                     // gamma = 4 -> -14
p3    = real(sd) - imag(sd)/tand(angd(sd - z3) - phi3);     // -20.0465
kc3   = 1/abs(horner(gcanc, sd)*(sd - z3)/(sd - p3));       // 6.4786e7
mprintf("phi=%.2f  z=%.2f  p=%.4f  kc=%.4e\n", phi3, z3, p3, kc3);

c_b    = syslin('c', kc3*(s + 1)*(s^2 + 2*s + 5)*(s - z3), ..
                      s*(s + 81)^2*(s - p3));
FTMF_b = (c_b*g) /. 1;

L_b     = syslin('c', kc3*(s - z3), s*(s + 81)^2*(s - p3)*(s + 15)^3);
polos_b = roots(s*(s + 81)^2*(s - p3)*(s + 15)^3 + kc3*(s - z3));
disp(polos_b);
disp(horner(FTMF_b, 0));                        // ganho DC = 1

scf(3); clf();
evans(L_b/kc3, 1.5*kc3);
xtitle('LR com o compensador final', 'Re(s)', 'Im(s)');
xs2png(3, 'Imagens/fig_ex4b_lr_final.png');

t   = 0:0.001:3;
y_b = csim('step', t, FTMF_b);
wn  = sig/cosd(bet);
s2  = syslin('c', wn^2, s^2 + 2*cosd(bet)*wn*s + wn^2);
y_2 = csim('step', t, s2);
idx = find(abs(y_b - 1) > 0.02);
mprintf("Mp real = %.2f%%   Mp 2a ordem = %.2f%%   ts(2%%) = %.3f s\n", ..
        (max(y_b) - 1)*100, (max(y_2) - 1)*100, t(max(idx)));

scf(4); clf();
plot(t, y_b, 'b'); plot(t, y_2, 'k--'); xgrid();
legend(['malha fechada (compensador final)'; '2a ordem pura'], 4);
xtitle('Item (b) - resposta ao degrau', 't [s]', 'y(t)');
xs2png(4, 'Imagens/fig_ex4b_resposta.png');

scf(5); clf();
plot_omega(3, 40, polos_b, [z3], -90, 12);
plot([-sig -sig], [imag(sd) -imag(sd)], 'k*');
xtitle('Plano s - regiao Omega (item b)', 'Re', 'Im');
xs2png(5, 'Imagens/fig_ex4b_omega.png');

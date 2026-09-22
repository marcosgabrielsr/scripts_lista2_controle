// Lista 2 - Exercicio 1 (b)
// Definindo Variáveis do Sistema
m = 100;
k = 400;
b = 240;

// Definindo Variaveis do Espaço de Estados
A = [0 -k/m; 1 -b/m];
B = [k/m; b/m];
C = [0 1];
sys = syslin('c',A,B,C);

// Executando simulacao
t = 0 : .01 : 10;
[y,x] = csim('step',t,sys);

// Visualizando variaveis de estado x1 e x2
plot(t,x(1,:),'b',t,x(2,:),'r');xgrid

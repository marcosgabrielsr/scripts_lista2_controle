// Lista 2 - Exercicio 1 (c)
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

// Plotando a Velocidade (em função da entrada degrau e das variaveis de estado)
// velocidade = y_dot = x2_dot = x1 - b/m*x2 + b/m*u
u = ones(1,length(t));
v = x(1,:) - (b/m)*x(2,:) + (b/m)*u;
plot(t,v,'g');xgrid

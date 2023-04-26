%%
%N = 5
%SR1

clear
clc

%Definiendo variables simbólivas
N=5;
x=sym('x',[1 N+1]);
    

VAR = x(1);
for i = 2:N
    VAR = [VAR x(i)];
end
VAR;

%Definiendo la función
f=0;
for i = 1:N-1
    f = f + (x(i) - 2*x(i+1)^2)^(2);
end
f;

%Definiendo punto inicial
x0 = [1 1 1 1 1];

%Evaluando en el método de Newton
tic
[TAB Xk] = SR1(f, VAR, x0)
toc

%%
%N = 10
%SR1

clear
clc

%Definiendo variables simbólivas
N=10;
x=sym('x',[1 N+1]);


VAR = x(1);
for i = 2:N
    VAR = [VAR x(i)];
end
VAR;

%Definiendo la función
f=0;
for i = 1:N-1
    f = f + (x(i) - 2*x(i+1)^2)^(2);
end
f;

%Definiendo el punto inicial
x0=zeros(1,N)+1;

%Evaluando en el método de Newton
tic
[TAB Xk] = SR1(f, VAR, x0)
toc
 

%%
%N = 100
%SR1

clear
clc

%Definiendo variables simbólivas
N=100;
x=sym('x',[1 N+1]);


VAR = x(1);
for i = 2:N
    VAR = [VAR x(i)];
end
VAR;

%Definiendo la función
f=0;
for i = 1:N-1
    f = f + (x(i) - 2*x(i+1)^2)^(2);
end
f;

%Definiendo punto inicial
x0=zeros(1,N)+1;

%Evaluando en el método de Newton
tic
[TAB Xk] = SR1(f, VAR, x0)
toc


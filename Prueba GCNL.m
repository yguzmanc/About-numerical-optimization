%% Taller 2
%Punto 2
%Inciso a)
%GCNL

clear
clc

syms x1 x2;
f(x1,x2)=exp(x1+3*x2-0.1)+exp(x1-3*x2-0.1)+exp(-x1-0.1);
VAR=[x1 x2];
x0 = [0.1 0.1];

tic
[TAB] = GCNL(f, VAR, x0)
toc

%%

%Inciso b)
%GCNL

clear
clc

syms x1 x2 x3;
f(x1,x2,x3)=x1^2 - 2*x1*x2 + 2*x2^2 - 2*x2 + x3^2 - x1*x3;
VAR=[x1 x2 x3];
x0 = [0.01 0.01 0.01];

tic
[TAB] = GCNL(f, VAR, x0)
toc

%%
%Inciso c)
%GCNL

clear
clc

syms x1 x2 x3 x4;
f(x1,x2,x3,x4)=x1^2 + x2^2 +2*x3^2 + x4^2 - 5*x1 - 5*x2 -21*x3 + 7*x4;
VAR=[x1 x2 x3 x4];
x0 = [0.1 0.1 0.1 0.1];

tic
[TAB] = GCNL(f, VAR, x0)
toc


%%
%Inciso d) N = 10
%GCNL

clear
clc

N=10;
x=sym('x',[1 N+1]);


VAR = x(1);
for i = 2:N+1
    VAR = [VAR x(i)];
end
VAR;


f=0;
for i = 1:N
    f = f + (x(i)^2)^(x(i+1)^(2) + 1) + (x(i+1)^2)^(x(i)^(2) + 1);
end
f;

x0=zeros(1,N+1)+0.0001;

tic
[TAB] = GCNL(f, VAR, x0)
toc

%%
%Inciso d) 
% N = 100
%GCNL

clear
clc

N=100;
x=sym('x',[1 N+1]);


VAR = x(1);
for i = 2:N+1
    VAR = [VAR x(i)];
end
VAR;


f=0;
for i = 1:N
    f = f + (x(i)^2)^(x(i+1)^(2) + 1) + (x(i+1)^2)^(x(i)^(2) + 1);
end
f;

x0=zeros(1,N+1)+0.0001;

tic
[TAB] = GCNL(f, VAR, x0)
toc
 

%%
%inciso e)
%GCNL

clear
clc

Q = [5 2 1
    2 7 3
    1 3 9];

c = [-9
    0
    -8];

syms x1 x2 x3;

x = [x1
    x2
    x3];

f=(1/2)*x.'*Q*x-c.'*x;
VAR=[x1 x2 x3];
x0 = [0.1 0.1 0.1];

tic
[TAB] = GCNL(f, VAR, x0)
toc


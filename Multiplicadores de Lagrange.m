%Multiplicadores de Lagrange
%%
%%%%Ninguna activa
clc
clear

syms x y

f(x,y) = y;
g(x,y) = 1 - x^2-y^2;
h(x,y) = x - y^2;
k(x,y) = x + y;

L(x,y,l2,l3) = f(x,y) + l2*h(x,y) + l3*k(x,y);

gra = jacobian (L(x,y,l2,l3),[x,y,l2,l3]);

[px,py,pl1,pl2] = solve(gra,x,y,l2,l3)

% hes(x,y,l1,l2,l3) = hessian(L(x,y,l1,l2,l3),[x,y,l1,l2,l3]);

% det(hes(px(1),py(1),pl(1)));
% 
% det(hes(px(2),py(2),pl(2)));
% 
% det(hes(px(3),py(3),pl(3)));

%%

%%% La primera activa

clc
clear

syms x y l1

f(x,y) = y;
g(x,y) = 1 - x^2-y^2;
h(x,y) = x - y^2;
k(x,y) = x + y;

L(x,y,l1) = f(x,y) + l1*g(x,y);

gra = jacobian (L(x,y,l1),[x,y,l1]);

[px,py,pl1] = solve(gra,x,y,l1)

hes(x,y,l1) = hessian(L(x,y,l1),[x,y,l1])

det(hes(px(1),py(1),pl1(1)))

det(hes(px(2),py(2),pl1(2)))

%det(hes(px(3),py(3),pl(3)));

%%

%%% La segunda activa

clc
clear

syms x y l2

f(x,y) = y;
g(x,y) = 1 - x^2-y^2;
h(x,y) = x - y^2;
k(x,y) = x + y;

L(x,y,l2) = f(x,y) + l2*h(x,y);

gra = jacobian (L(x,y,l2),[x,y,l2]);

[px,py,pl2] = solve(gra,x,y,l2)

hes(x,y,l2) = hessian(L(x,y,l2),[x,y,l2])

det(hes(px(1),py(1),pl2(1)))

det(hes(px(2),py(2),pl2(2)))

%%
%%%Segundo y tercero activo
clc
clear

syms x y l3 l2

f(x,y) = y;
%g(x,y) = 1 - x^2-y^2;
h(x,y) = x - y^2;
k(x,y) = x + y;

L(x,y,l2,l3) = f(x,y) + l2*h(x,y) + l3*k(x,y);

gra = jacobian (L(x,y,l2,l3),[x,y,l2,l3]);

[px,py,pl2,pl3] = solve(gra,x,y,l2,l3)

hes(x,y,l2,l3) = hessian(L(x,y,l2,l3),[x,y,l2,l3])

det(hes(px(1),py(1),pl2(1),pl3(1)))
det(hes(px(2),py(2),pl2(2),pl3(2)))

%%

%%%La primera y segunda activa
clc
clear

syms x y l1 l2

f(x,y) = y;
g(x,y) = 1 - x^2-y^2;
h(x,y) = x - y^2;
k(x,y) = x + y;

L(x,y,l1,l2) = f(x,y) + l1*g(x,y) + l2*h(x,y);

gra = jacobian (L(x,y,l1,l2),[x,y,l1,l2]);

[px,py,pl1,pl2] = solve(gra,x,y,l1,l2)

hes(x,y,l1,l2) = hessian(L(x,y,l1,l2),[x,y,l1,l2])

det(hes(px(1),py(1),pl1(1),pl2(1)))
det(hes(px(2),py(2),pl1(2),pl2(2)))
det(hes(px(3),py(3),pl1(3),pl2(3)))
det(hes(px(4),py(4),pl1(4),pl2(4)))

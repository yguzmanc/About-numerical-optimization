%%%%%%%%%%%%%%%%Método del Steepest Descent%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Definiendo la función, su gradiente y demás datos iniciales
clear
clc

%Instrucción para escribir la función y se almacena como una cadena de
%carácteres.
f = input('Escriba la función','s');

%Convertir la funcion a simbólica
f = str2sym(f);

%disp('Gradiente:')
%Determinar el vector de variables
Nv = symvar(f);

%Calcular el gradiente
grad = gradient(f,Nv);

%Clacular el Hessiano
hess = jacobian(grad,Nv);

%Mostrar componentes del gradiente
%disp(grad);

% ¿Se Minimiza o maximiza?
M = input('Si se busca el mínimo de la función ingrese el número 0, si se busca el máximo ingrese el 1');

%¿Cuál es el punto inicial x0?
x0 = input('Punto inicial (ingrese un vector fila)');

%¿Cuál es la tolerancia?
tol = input('Defina la tolerancia');

%Máximo número de iteraciones
IterMax = input('Establezca el número de iteraciones máximas en caso de no alcanzar la tolerancia.');

%%
%Inicio del algoritmo

Egrad = subs(grad,Nv,x0);

Ehess = subs(hess,Nv,x0);

if M == 1
    Egrad=Egrad;
else
    Egrad=-Egrad;
end

Iter = 1;

while norm(Egrad)>=tol && Iter<=IterMax
    [alpha]=Newton(grad,hess,x0)

    D(Iter,:) = {Iter,x0,Egrad,alpha}

    x0 = x0+alpha*Egrad;

    Egrad = subs(grad,Nv,x0);

    Ehess = subs(hess,Nv,x0);

    if M == 1
    Egrad=Egrad;
    else
    Egrad=-Egrad;
    end
    Iter = Iter + 1;

    D(Iter,:) = {Iter,x0,Egrad,alpha}
end

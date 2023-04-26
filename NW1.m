function [TAB, Xk] = NW1(f, VAR, x0)

% Calculo del gradiente y el jacobiano:
VAR=VAR.';
G = (gradient(f,VAR));
Geva = [1e10,1e10];

J = jacobian(G,VAR);

e = 10^(-2);
i=0;
Xk = x0.';
x_datos=[];
Px_datos=[];
norm_datos = [];

% Condiciones de optimizacion:
while norm(Geva) > e
    Geva = subs(G,VAR,Xk); %Gradiente actualizado
    Jeva = subs(J,VAR,Xk);%Jacobiano actualizado
    S = inv(double(Jeva)); %Nueva búsqueda de direccion
    Pk = -S*Geva;
    Xk = Xk+Pk;
    Xk = double(Xk);
    x_datos = [x_datos Xk];
    Px_datos = [Px_datos Pk];
    norm_datos = [norm_datos norm(Geva)];
    TAB = [x_datos.',Px_datos.',norm_datos.']
    i=i+1;


end

%Tabla de resultados:`
Iter = 1:i;


Iteraciones = Iter';


TAB = [Iteraciones, x_datos.',Px_datos.',norm_datos.'];
TAB = double(TAB);


if (norm(Geva) < e)
    fprintf('Mínimo obtenido satisfactoriamente...\n\n');
end

end

%***************** WOLFE CONDITION LINE SEARCH *****************%
% ---------------------------------------------------------------
% Función para generar alpha, que comple las condiciones de Wolfe
% ---------------------------------------------------------------
% Parámetros:
%   f: Función en forma simbólica
%   df: Gradiente de f en forma simbólica
%   V: Vector con las variables de f
%   xk: Vector con el valor de la iteración anterior
%   p: Dirección del descenso
% Retorno:
%   a: valor de alpha que cumple las condiciones de Wolfe
function a_star = WolfeConditionLineSearch(f, V, xk, p)
    a0 = 0;
    a_max = random(2,5);    % Valor aleatoreo mayor a cero
    a1 = random(0, a_max);
    c1 = 10^(-4);           % Valor sugerido para 'c1'
    c2 = random(c1,1);
    i = 1;
    % Definir la función 'phi' y su derivada:
    syms a;
    phi = eval(subs(f, V, xk + a*p));
    dphi = diff(phi, a);
    while abs(a1 - a0) >= 10^(-6)
        phi_a0 = eval(subs(phi, a, a0));
        phi_a1 = eval(subs(phi, a, a1));
        phi_0  = eval(subs(phi, a, 0));
        dphi_0 = eval(subs(dphi, a, 0));
        if (phi_a1 > phi_0 + c1 * a1 * dphi_0) || (phi_a1 >= phi_a0 && i > 1)
            a_star = zoom(phi, dphi, a, a0, a1, c1, c2);
            break;
        end
        dphi_a1 = eval(subs(dphi, a, a1));
        if abs(dphi_a1) <= -c2 * dphi_0
            a_star = a1;
            break;
        end
        if dphi_a1 >= 0
            a_star = zoom(phi, dphi, a, a1, a0, c1, c2);
            break;
        end
        a0 = a1;
        a1 = random(a1, a_max);
        i = i + 1;
        a_star = a1;
    end
return 

% ------------------------ Fin de WolfeConditionLineSearch

function a_star = zoom(phi, dphi, a, a0, a1, c1, c2)
    aj = 0; aprev = 1;
    while abs(aj - aprev) >= 10^(-6)
        phi_a0 = eval(subs(phi, a, a0));
        phi_a1 = eval(subs(phi, a, a1));
        phi_0 = eval(subs(phi, a, 0));
        dphi_0 = eval(subs(dphi, a, 0));
        % Interpolación
        aprev = aj;
        % Interpolación cúbica:
        dphi_a0 = eval(subs(dphi, a, a0));
        dphi_a1 = eval(subs(dphi, a, a1));
        if (a0 < a1)
            aj = InterpolacionCubica(a0, a1, phi_a0, phi_a1, dphi_a0, dphi_a1);
        else
            aj = InterpolacionCubica(a1, a0, phi_a1, phi_a0, dphi_a1, dphi_a0);
        end
        %aj = (a0 + a1)/2;
        % -------------------------------------------------------
        phi_aj = eval(subs(phi, a, aj));
        if (phi_aj > phi_0 + c1 * aj * dphi_0) || (phi_aj >= phi_a0)
            a1 = aj;
        else
            dphi_aj = eval(subs(dphi, a, aj));
            if abs(dphi_aj) <= -c2*dphi_0
                a_star = aj;
                break;
            end
            if dphi_aj * (a1-a0) >= 0
                a1 = a0;
            end
            a0 = aj;
        end
    end
    a_star = aj;
return % -------------------------------------------- Fin de Zoom


% ---------------------------------------------------------------
% Función para generar un número aleatorio en el intervalo [a, b]
% ---------------------------------------------------------------
% Parámetros:
%   a: Extremo izquierdo del intervalo
%   b: Extremo derecho del intervalo
% Retorno:
%   x: Valor aleatareo, que cumple que a < x < b
function x = random(a, b)
    a = a;
    b = b;
    x = rand(1) * (b - a) + a;
return % ------------------------------------------ Fin de random


% Algoritmo que calcula 'alpha_j' con interpolacion cúbica
% @author: http://www.math.unl.edu/~tshores1/Public/OctaveFiles/OptimizeTools/cubicmin.m
function aj = InterpolacionCubica(a, b, phi_a, phi_b, dphi_a, dphi_b)
    % revisar página 59 del libro.
    h = b - a;
    phi_ab = (phi_b - phi_a)/h;
    % Calculo los d´s (Page 59, PDF-page 80)
    d1 = dphi_a + dphi_b - 3*phi_ab;
    d2 = sqrt(d1*d1 - dphi_a*dphi_b);
    absc = [a b];
    ord = [phi_a, phi_b];
    if (isreal(d2)) % Ver si d2 no es imaginario
        phi_aab = (phi_ab-dphi_a)/h;
        phi_aabb = (dphi_b-2*phi_ab+dphi_a)/(h*h);
        if (abs(phi_aabb) < 100*eps)
            if (abs(phi_aab) < 100*eps)
                x = a;
            else
                x = a - dphi_a/(2*phi_aab);
            end
        else
            x = b - h*(dphi_b + d2 - d1)/(dphi_b - dphi_a + 2*d2);
        end
        x = max(min(x,b),a);
        absc = [a b x];
        ord = [phi_a,phi_b, phi_a + (x-a)*(dphi_a + (x-a)*(phi_aab + (x-b)*phi_aabb))];
    end
    aj = min(ord);
    aj = find(ord == aj);
    aj = absc(aj(1));
return % ----------------------------- Fin de InterpolacionCubica
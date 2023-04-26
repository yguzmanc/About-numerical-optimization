function[TAB] = GCNL(f,VAR,x0)
%Gradiente Conjugado No Lineal
%En particular "The Fletcher-Reeves method"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f0 = eval(subs(f,VAR,x0));
gradf = gradient(f);
gradf0 = eval(subs(gradf,VAR,x0));
gradfk=gradf0;
p0 = -gradf0';
pk = p0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estableciendo los parámetro
i = 1;
max_iter = 100;
norm_k = norm(gradf0);
ak = WolfeConditionLineSearch(f, VAR, x0, p0);
xk_plus1 = x0+ak*pk;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Iniciando las iteraciones
while (norm_k>=10^(-3) && i<=max_iter)
    gradfk_plus1 = eval(subs(gradf,VAR,xk_plus1));
    Bk_plus1 = (gradfk_plus1'*gradfk_plus1)/(gradfk'*gradfk);
    pk_plus1 = -gradfk_plus1+Bk_plus1*pk';
    pk_plus1 =  pk_plus1';
    ak= WolfeConditionLineSearch(f, VAR, xk_plus1, pk_plus1);
    xk_plus1 = xk_plus1+ak*pk_plus1;
    norm_k = norm(gradfk_plus1)
    %pk_plus1 Le quité esto a TAB
    TAB(i,:) = [i xk_plus1 norm_k]
    i=i+1;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
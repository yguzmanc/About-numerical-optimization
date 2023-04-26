function[TAB Xk] = SR1(f,V,x0)


%Initialize tolerance
Tol=10^(-3);
IteMax = 1000;

I = eye(length(V));     
df = jacobian(f,V);
df=df';

%Initialize xk, Hk
Xk = x0';
Hk = I;
grad_Xk = eval(subs(df,V,Xk'));
grad_Xk = grad_Xk;
e = norm(grad_Xk);
K=0;
OutputTable = [];

while (e > Tol) && (K<IteMax)
    % Determine alpha k
    p_k = -Hk*grad_Xk;
    alphaK = WolfeConditionLineSearch(f,V,Xk',p_k');
    % Update x_k and find new gradient
    x_new = Xk + alphaK*p_k;
    ngrd = eval(subs(df,V,x_new'));
   
    % Update inverse Hessian Hk+i
    Sk = x_new - Xk;
    Yk = ngrd - grad_Xk;
    pk_h = 1/(Yk'*Sk);
    coeff = Sk-Hk*Yk;
    Hk = Hk + (coeff*coeff')/(coeff'*Yk);
    Xk = x_new;
    grad_Xk = ngrd;
    e = norm(grad_Xk);
    K = K + 1;
    % Store the result in a table
    OutputTable(K,:) = [ K  Xk' p_k' e];
    formatSpec = 'Iteracion %d y error: %f\n';
    %str = sprintf(formatSpec,K,e);
    %disp(str);
end

    TAB = OutputTable;

end
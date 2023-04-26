function[TAB,Xk] = BFGS(f,VAR,x0)

Tol=10^(-3);

VAR = VAR.';

IteMax = 150;

IdenMatrix = eye(length(VAR));     
gradf = jacobian(f,VAR).';

Xk = x0.'
Hk = IdenMatrix;
grad_Xk = eval(subs(gradf,VAR,Xk));
e = norm(grad_Xk);
K=0;
OutputTable = [];

while (e > Tol) && (K<IteMax)
    % Determine alpha k
    p_k = -Hk*grad_Xk;
    alphaK = WolfeConditionLineSearch(f,VAR,Xk,p_k);
    % Update x_k and find new gradient
    x_new = Xk + alphaK*p_k;
    newGrad = eval(subs(gradf,VAR,x_new));
   
    % Update inverse Hessian Hk+i
    Sk = x_new - Xk;
    Yk = newGrad - grad_Xk;
    pk_h = 1/(Yk.'*Sk);
    Hk = (IdenMatrix-pk_h*(Sk*Yk'))*Hk*(IdenMatrix-pk_h*(Yk*Sk')) + pk_h*(Sk*Sk');
    Xk = x_new;
    grad_Xk = newGrad;
    e = norm(grad_Xk);
    K = K + 1;
    % Store the result in a table
    OutputTable(K,:) = [ K Xk' p_k' ];
    formatSpec = 'Iteracion %d y error: %f\n';
    str = sprintf(formatSpec,K,e);
    disp(str);
end

    disp('f converge a Xk en ');
    Xk;
    TAB = OutputTable(1:K,:);
end
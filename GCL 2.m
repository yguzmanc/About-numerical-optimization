function [TAB,Xk] = GCL(A,b,x0)
maxIterations = 200;
Tolerance = 10^(-6);
i = 0;
TAB = zeros(maxIterations, 1 + 2*length(b)+1) ;

rk = A*x0.'-b.';
Xk= x0.';
Pk = -rk;
while norm(rk)>Tolerance && (i<maxIterations)
    i = i+1;
    alphak = (rk.'*rk) / (Pk.'*A*Pk);   
    Xk = Xk+alphak*Pk;
    rk1 = rk + alphak*A*Pk;
    bk = (rk1.'*rk1) / (rk.'*rk);
    TAB(i,:) = [i Xk.'  Pk.'  norm(rk1)];
    rk = rk1;
    Pk = -rk1+bk*Pk;
    disp(Pk');
end

for j = i+1:maxIterations
    TAB(i+1,:) = [];
end

end
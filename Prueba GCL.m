clear
clc

A = [1 0 0
    0 2 0
    0 0 3];

b = [1 1 1];

x0 = [0 0 0];

[TAB, Xk] = GCL(A,b,x0)
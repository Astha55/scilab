               // Preliminay exercise 1
// Using built in functions

clc
clear
N=7
V=2
C=3

// Second order derivative matrix
A=-2*eye(N-2,N-2)+diag(ones(N-3,1),1)+diag(ones(N-3,1),-1)
disp("SD matrix",A)

// kinetic energy matrix
B=C*A
disp("KE",B)

//Potential Energy matrix
v=V*ones(1,N-2)
P=diag(v)
disp("PE",P)

//Hamiltonian matrix
H=B+P
disp("Hamiltonian matrix", H)


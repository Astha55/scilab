             // Preliminary exercise 1
// Using loops

clc
clear
n=7
v=2
c=3

// Second derivative matrix
for i=1:n-2
    for j=1:n-2
        if abs(i-j)==1
            A(i,j)=1
        elseif i==j
            A(i,j)=-2
        end
    end
end
disp("SD matrix",A)

// Kinetic Energy matrix
B=c*A
disp("KE",B)

//Potential energy matrix
for i=1:n-2
    P(i,i)=v
end
disp("PE",P)

// Hamiltonian Matrix
H=B+P
disp("Hamiltonian Matrix",H)


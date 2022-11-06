// Date: 13-09-2022
// Program 3: Sensitivity of solutions to the choice of parameters.
// Varying the width of the potential well.
clc
clear

Lscale = 1D-9;
e=1.6*10^-19
h=1.054*10^-34
m=9.109*10^-31
n=1001
U=-100

t(1)=0.05
t(2)=0.2
t(3)=0.4
r(1)=0.1
r(2)=0.2
r(3)=0.4

for j=1:3
        disp("For width, t = ",t(j))
        x1=-r(j)*Lscale
        xn=+r(j)*Lscale
        x=linspace(x1,xn,n)
        del_x=(xn-x1)/(n-1)
        //disp("The step size is = ",del_x)
        C=(1/e)*((-h^2)/(2*m*(del_x^2)))
        //disp("The C is =",C)
        
        A=zeros(n-2,n-2)
        // By Using iterative loops
        for i=1:n-2
            if i==1
                A(i,i+1)=1
            elseif i==n-2
                A(i,i-1) =1
            end
            A(i,i)=-2
        end
        for i=2:n-3
                A(i,i-1)=1
                A(i,i+1)=1
        end
        K=C*A
        
        // defining the matrix for variable potential
        function y=fun(x)
            l=length(x)
            for i=1:del_x:l
                if x(i)>=-t(j)*Lscale && x(i)<=t(j)*Lscale then
                    y(i)=U
                else
                    y(i)=0
                end
            end
        endfunction
        for i=1:n
                v(i)=fun(x(i))
        end
        V= diag(v(2:n-1))
        H=K+V
        [Evec,Eval]=spec(H)
        
        counter=0
        for i=1:n-2
            if Eval(i,i)<0
                counter=counter+1
            end
        end
        disp("The number of bound states is = ",counter)
        
        disp("The Energy Eigen values are = ")
        for i=1:counter
            P(i,1)=Eval(i,i)
        end
        disp(P)
end


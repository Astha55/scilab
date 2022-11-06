// Date: 13 Sept, 2022
// Program 3: Sensitivity of solutions to the choice of parameters.
// Variation in the Value of N.

clc
clear

Lscale= 1D-9;
x1= -0.2*Lscale      //min
xn= +0.2*Lscale      //max
e=1.6*10^-19
h=1.054*10^-34
m=9.109*10^-31

n(1)=11       //comparing w small

n(2)=10
n(3)=100
n(4)=501
n(5)=1001
n(6)=2001

for j=1:6         //Starting loop
        disp(" For n= ",n(j))
        x=linspace(x1,xn,n(j))
        del_x=(xn-x1)/(n(j)-1)
        //disp("The step size is = ",del_x)
        
        
        C=(1/e)*((-h^2)/(2*m*(del_x^2)))
        //disp("The C is =",C)
        
        A=zeros(n(j)-2,n(j)-2)
        // By Using iterative loops
        for i=1:n(j)-2
            if i==1
                A(i,i+1)=1
            elseif i==n(j)-2
                A(i,i-1) =1
            end
            A(i,i)=-2
        end
        for i=2:n(j)-3
                A(i,i-1)=1
                A(i,i+1)=1
        end
        K=C*A
        
        // defining the matrix for variable potential
        function y=fun(x)
            l=length(x)
            for i=1:del_x:l
                if x(i)>= -0.1*Lscale && x(i)<=0.1*Lscale then
                    y(i)= -100
                else
                    y(i)=0
                end
            end
        endfunction
        for i=1:n(j)
                v(i)=fun(x(i))
        end
        V= diag(v(2:n(j)-1))
        
        H=K+V
        [Evec,Eval]=spec(H)
        
        counter=0
        for i=1:n(j)-2
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

// Date: 15-08-2022
// Program 1: Eigenvalues and Eigenvectors.

clc
clear
clf

Lscale = 1D-9;
x1=-0.2*Lscale
xn=+0.2*Lscale
n=1001
U=100
x=linspace(x1,xn,n)
del_x=(xn-x1)/(n-1)
//disp("The step size is = ",del_x)

e=1.6*10^-19
h=1.054*10^-34
m=9.109*10^-31
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
        if x(i)>=-0.1*Lscale && x(i)<=0.1*Lscale then
            y(i)=-U
        else
            y(i)=0
        end
    end
endfunction
for i=1:n
        v(i)=fun(x(i))
end
V= diag(v(2:n-1))
axis=gca()
axis.data_bounds=[x1,-(U+15);xn, 15]
plot(x',v)
title("Plot of (U(x) Vs x) ",'fontsize',5,'color','r')


H=K+V
[Evec,Eval]=spec(H)

counter=0
for i=1:n-2
    if Eval(i,i)<0
        counter=counter+1
    end
end
disp("The number of bound states is = ",counter)

disp("The Bound state Eigen values are = ")
for i=1:counter
    S(i,1)=Eval(i,i)
end
disp(S)

// Plot of U(x) Versus x - Energy level diagram
e1=S(1,1)*ones(n,1)
plot(x',e1,'-dr')
e2=S(2,1)*ones(n,1)
plot(x',e2,'-^g')
e3=S(3,1)*ones(n,1)
plot(x',e3,'-sc')
e4=S(4,1)*ones(n,1)
plot(x',e4,'-+k')

xnew=x(2:n-1)
xnewt=xnew'

for i=1:counter
    c=(Evec(:,i).*Evec(:,i))
    P(i)=inttrap(xnew,c)
    Q(:,i)=(1/sqrt(P(i)))*Evec(:,i)
end

show_window(1)
for i=1:counter
    subplot(2,2,i)
    plot(xnew',Q(:,i))
    title("Plot - Normalized Eigenfunctions Vs x for the bound state - "+string(i),'fontsize',5,'color','g')
    xlabel('------ x axis -------->','Fontsize',5,'color','orange')
    ylabel('------ y axis -------->','Fontsize',5,'color','orange')
    xgrid
end

show_window(2)
for i=1:counter
    c(:,i)=(Evec(:,i).*Evec(:,i))
    subplot(2,2,i)
    plot(xnew',c(:,i))
    title("Plot -  |psi(x)^2| Vs x for the bound state - "+string(i),'fontsize',5,'color','g')
    xlabel('------ x axis -------->','Fontsize',4,'color','brown')
    ylabel('------ y axis -------->','Fontsize',4,'color','brown')
    xgrid
end

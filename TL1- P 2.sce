// Date: 15-08-2022
// Program 2: Expectation values and other physical info.

clc
clear

Lscale = 1D-9;
x1=-0.2*Lscale
xn=+0.2*Lscale
n=1001
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
            y(i)=-100
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
//disp("The number of bound states is = ",counter)

//disp("The first four Eigen values are = ")
for i=1:counter
    ng(i,1)=Eval(i,i)
end
//disp(ng)

xnew=x(2:n-1)
xnewt=xnew'

for i=1:counter
    c=(Evec(:,i).*Evec(:,i))
    P(i)=inttrap(xnew,c)
    Q(:,i)=(1/sqrt(P(i)))*Evec(:,i)
end

// Expectation value of x
for i=1:counter
    r=(Q(:,i).*Q(:,i))
    nv3(:,i)=(xnewt.*r)
    Ex_x(i)=inttrap(xnew,nv3(:,i))
    end
disp('Expectation value of x for all bound states are = ',Ex_x)

// Expectation value of x^2
for i=1:counter
    Xnew=(xnew.*xnew)
    s=Q(:,i).*Q(:,i)
    nv4(:,i)=((Xnew').*s)
    Ex_xx(i)=inttrap(xnew,nv4(:,i))
    end
disp('Expectation value of x^2 for all bound states are = ',Ex_xx)

Vnew=v(2:n-1)
// Expectation value of V
for i=1:counter
    t=Q(:,i).*Q(:,i)
    nv5(:,i)=((Vnew).*t)
    Ex_V(i)=inttrap(xnew,nv5(:,i))
    end
disp('Expectation value of V for all bound states are = ',Ex_V)

// Expectation value of p
for i=1:n-3
    for j=1:counter
        mid_s(i,j)=(Q(i+1,j)+Q(i,j))/2
        diff_s(i,j)=(Q(i+1,j)-Q(i,j))/del_x
        Xmid(i)=(xnewt(i)+xnewt(i+1))/2
    end    
end
for j=1:counter
    y2(:,j)=mid_s(:,j).*diff_s(:,j)
    Ex_p(j)=(-(%i*h))*(inttrap(Xmid,y2(:,j)))
    end
disp('Expectation value of p for all bound states are = ',Ex_p)

// Expectation value of p^2
for i=1:n-4
    for j=1:counter
        mid_s2(i,j)=(mid_s(i+1,j)+mid_s(i,j))/2
        diff_s2(i,j)=(Q(i+2,j)-2.*Q(i+1,j)+Q(i,j))/(del_x^2)
        Xmid2(i)=(Xmid(i)+Xmid(i+1))/2
    end    
end
for j=1:counter
    y22(:,j)=mid_s2(:,j).*diff_s2(:,j)
    Ex_p2(j)=(-(h^2))*(inttrap(Xmid2,y22(:,j)))
    
end
disp('Expectation value of p^2 for all bound states are = ',Ex_p2)
    
for i=1:counter
    Ex_ke(i)=(Ex_p2(i)/(2*m*e))
    Te(i)=Ex_V(i)+Ex_ke(i)
end
disp('Expectation value of K E for all bound states are = ',Ex_ke)
disp('Expectation value of Te for all bound states are = ',Te)
    
    //Uncertainty Product       
Vx=Ex_xx-(Ex_x)^2
Vp=Ex_p2-(Ex_p2)^2
sigma_x=sqrt(Vx)
sigma_p=sqrt(Vp)
u=(sigma_p.*sigma_x)/(h/2)
disp("Uncertainity product = ",u)

                 //Preliminary Exercise 3
// Energy eigen values for Finite Square Well Potential
clc
clear
clf

h=6.63D-34         //value of planks constant
m=9.1D-31        //mass of an electron
e=1.6D-19        //charge on an electron
xmin=-0.2
xmax=0.2
N=1001
vt=0.1
vo=-100
x=linspace(xmin,xmax,N)
dx=(x(2)-x(1))*1D-9
c=-((h/(2*%pi))^2)/(2*m*(dx^2)*e)

a=ones(1,N-2)
b=ones(1,N-3)
V=zeros(1,N)
P=diag(-2*a)
Q=diag(1*b,1)
R=diag(1*b,-1)
T=P+Q+R

for i=1:N              //using loops
    if abs(x(i))<=vt then
        V(i)=vo
    end
end

v=diag(V(2:N-1))

M=c*T+v

xnew=x(2:N-1)

plot(xnew,V(2:N-1))

[a1,a2]=spec(M)
z=spec(M)


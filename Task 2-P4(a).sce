//Date: 28 sept,2022
//Aim: Conventional Hydrogenic Atoms
//Replacement of nuclear charge.
clc
clear
clf

N=1001
H=6.626D-34
m=9.1D-31
x0=0.1D-14
xm=1D-9
h=(xm-x0)/(N-1)
x=x0:h:xm
e=1.6D-19
E0=9D9                 //Coulomb Constant
C=-((H/(2*%pi))^2)/(2*m*h*h*e)
z(1)=1
z(2)=2
z(3)=3
for l=1:3
    disp('for z=',l)
    v=-((E0*e*(z(l)^2)./x))
    P=diag(((-2*C)*ones(N-2,1)))+diag(((C)*ones(N-3,1),-1))+diag(((C)*ones(N-3,1),1))
    Q=diag(v((2:N-1)))
    R=P+Q
    [E,J]=spec(R)
    x_new=x(2:N-1)
    s=0;
    for i=1:N-2
        if J(i,i)<0
            s=s+1
        end
    end
    disp("No. of bound states : ",s)
    disp('<--Energy eigenvalue for Bound States-->')
    for i=1:s
        disp(J(i,i))
    end
    disp('Bohr theo. energy')
    for i=1:s
        B_En=-(13.6*((z(l))^2))/((i)^2)
        disp(B_En)
    end
    for n=1:s
        Y=E(:,n).*E(:,n)
        [a,b]=max(Y)
        Bohr_R=((n^2)*((H/(2*%pi))^2)/((9*1D9)*m*e^2))/z(l)
        disp('Bohr radius',Bohr_R)
        disp('r_mp',x(b))
    end
end

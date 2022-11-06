//Date: 28 sept,2022
//Aim: Muonic,Pionic and kaonic Hydrogen.
clc
clear
clf

N=1001
H=6.626D-34
m=9.1D-31
x0=0.1D-14
xm=4D-12
e=1.6D-19
h=(xm-x0)/(N-1)
m_u=[207,273,966]
m_VH=m_u*m
for k=1:3
    if(k==1)
        disp('For muonic')
    elseif(k==2)
        disp('for pionic')
    else
        disp('for kaonic')
    end
    Bohr_R=((H/(2*%pi))^2)/((9*1D9)*m_VH(k)*e^2)
    x=x0:h:xm
    E0=9D9                 //Coulomb Constant
    C=-((H/(2*%pi))^2)/(2*m_VH(k)*h*h*e)
    v=-((E0*e)./x)
    P=diag(((-2*C)*ones(N-2,1)))+diag(((C)*ones(N-3,1),-1))+diag(((C)*ones(N-3,1),1))
    Q=diag(v(2:N-1))
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
    disp('Bohr theo. energy')
    for i=1:s
        disp(J(i,i))
        B_En=-13.6/(i^2)
        disp(B_En)
    end
    for n=1:s
    Y=E(:,n).*E(:,n)
    [a,b]=max(Y)
    B_R=(n^2)*((H/(2*%pi))^2)/((9*1D9)*m_VH(k)*e^2)
    disp('Bohr radius',Bohr_R)
    disp('r_mp',x(b))
    end
    //Normalization loop
    for i=1:s
        Y=E(:,i).*E(:,i)
        I(i)=inttrap(x_new,Y)
        N_E(:,i)=(1/sqrt(I(i))).*E(:,i)
        r=x_new'
        U(:,i)=N_E(:,1)./r
    end
    for i=1:s
        show_window(1)
        subplot(2,2,i)
        title('<---Normalised Eigenvectors plot-->')
        plot(x_new',N_E(:,i))
    end
end

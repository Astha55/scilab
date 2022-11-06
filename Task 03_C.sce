// P3: Sensitivity to choice of Parameters.

clc
clear
clf

N=1001
H=6.626D-34                //Planck's Constant
m=9.1D-31                 //mass of electron
e=1.6D-19              //Charge on electron
x0=0.1D-14
xm(1)=1D-7
xm(2)=1D-9
xm(3)=1D-11
a=[0.3D-9,0.5D-9,0.7D-9]
E0= 9D9                 //Coulomb Constant

for s=1:3
    disp('<-----for rmax = ',s)
    h=(xm(2)-x0)/(N-1)        // For del x
    x=x0:h:(xm(s))
    C=-((H/(2*%pi))^2)/(2*m*h*h*e)
    for k=1:3
        v(k,:)=-((E0*e.*exp(-x/a(k)))./x)

        if(k==1)
            plot(x(3:N),v(3:N),':r')
            title('for rmax=1D-7')
        elseif(k==2)
            show_window(1)
            plot(x(3:N),v(3:N),':r')
            title('for rmax=1D-9')
        else
            show_window(2)
            plot(x(3:N),v(3:N),':r')
            title('for rmax=1D-11')
        end

    end
end


//Date: 14 sept, 2022
//Aim: Coulomb Potential for Hydrogen Atom and Hydrogenic Systems.

clc
clear
clf

N=1001
rmin=1*(1e-15)
rmax=1*(1e-9)
 
r=linspace(rmin,rmax,N)
del_r=(r(2)-r(1))
r_new=(r(2:N-1))'

a=0.529D-10                  //bohr's radius
Z=1                     //atomic no. of hydrogen atom.
ec=9e9              //electrostatic constant = 1/(4*%pi*e0)
h=6.63e-34    
m=9.109e-31    
e=1.6e-19     
C=-((h/(2*%pi))^2)/(2*m*(del_r^2)*e)  //value of the constant
h_cut=(h/(2*%pi))
//disp(h_cut)
//Second Derivative & KINETIC ENERGY 

C1=-2
D=zeros(N-2,N-2)

for i=1:N-2
    D(i,i)=C1
    if i==1
        D(i,i+1)=1
    elseif i==(N-2)
        D(i,i-1)=1
    else
        D(i,i+1)=1
        D(i,i-1)=1 
    end
end

KE=C*D         //Diagonal Matrix

//POTENTIAL ENERGY
V0=-((Z*e*ec)./r)        //Potential value in electron volt.
V=diag(V0(2:N-1))
H=KE+V         
//disp("spec(H)---------------",spec(H))
[eigfn,eigval]=spec(H)

//P1: Eigenvalues & Eigenvectors

show_window(0)
plot(r(3:N)',V0(3:N)')
title("Plot - Potential v/s r ","fontsize",5,"color","Purple")
xlabel('--------x-axis------------>''Fontsize',4,'color','Purple')
ylabel('--------y-axis------------>''Fontsize',4,'color','Purple')

///For calculating the no. of particle in the bound state of finite well we use counter:

counter=0
for i=1:N-2
    if eigval(i,i)<0
        counter=counter+1
    end
end
disp("COUNTER",counter)
disp("The first three eigen value are:")

for i=1:counter
    disp(eigval(i,i))
end

for i=1:counter
    nv=eigfn(:,i)
    norml(i)=(inttrap(r_new,(nv.^2)))     //Normalization of the bound state
    WF(:,i) = eigfn(:,i)./sqrt(norml(i))
    R(:,i)=((WF(:,i))./r_new)
end

for i=1:counter
    show_window(1)
    subplot(2,2,i) 
    plot(r_new,WF(:,i))             //xnew v/s wavefunction

    title("Plot - Eigenfunction v/s r for the bound state.-."+string(i),"fontsize",5,"color","green")
    
    xlabel('--------x-axis------------>''Fontsize',4,'color','purple')
    ylabel('--------y-axis------------>''Fontsize',4,'color','purple')
    g=gca()

    g.x_ticks=tlist(["ticks","locations","labels"],[5.010D-13,5.000D-10,9.995D-10],["5.010D-13","5.000D-10","9.995D-10"]);
    xgrid
end

for i=1:counter
    show_window(2)
    subplot(2,2,i) 
    plot(r_new,WF(:,i).^2)              //Probability Density
    
    title("Plot - Probability Density v/s r for the bound state.-."+string(i),"fontsize",5,"color","M")
    
    xlabel('--------x-axis------------>''Fontsize',4,'color','purple')
    ylabel('--------y-axis------------>''Fontsize',4,'color','purple')
    g=gca()

    g.x_ticks=tlist(["ticks","locations","labels"],[5.010D-13,5.000D-10,9.995D-10],["5.010D-13","5.000D-10","9.995D-10"])
    xgrid
end

///R(r) theoritical value:
Rth(:,1)=2*(a^(-3/2))*exp(-r_new./a)
Rth(:,2)=(-1/sqrt(2))*(a^(-3/2))*(1-(1/2*(r_new./a))).*exp(-r_new./(2*a))
Rth(:,3)=(-2/3*sqrt(3))*(a^(-3/2))*(1-(2/3*(r_new./a))+(2/27*((r_new./a)^2))).*exp(-r_new./(3*a))

for i=1:counter
    show_window(3)
    subplot(2,2,i) 
    plot(r_new,R(:,i))            //Another plot
    plot(r_new,Rth(:,i),"--r")

    title("Plot - R v/s r for the bound state.-."+string(i),"fontsize",5,"color","blue")
    xlabel('--------x-axis------------>''Fontsize',4,'color','black')
    ylabel('--------y-axis------------>''Fontsize',4,'color','black')
    g=gca()
    
    g.x_ticks=tlist(["ticks","locations","labels"],[5.010D-13,5.000D-10,9.995D-10],["5.010D-13","5.000D-10","9.995D-10"])
    xgrid
end

disp("The energy from Bohrs spectrum are:")
for i=1:counter
    E=-13.6/(i^2)
    disp(E)
end

///Most Probable value of r
for i=1:counter
    [p,q]=max(WF(:,i))
    r_mpv=r(q)
    disp("Most probable value of r:",r_mpv)
    r_n=(i^2).*a         //bohrs radius using formula
    disp(r_n)
end

///Verification of Orthonormal set:
for i=1:counter
    for j=1:counter
        if (i<>j)
            Oth(i)=inttrap(r_new,(WF(:,i).*WF(:,j)))
        end
    end
end
disp("Orthonormal set")
disp(Oth)






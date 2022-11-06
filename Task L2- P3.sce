//Date; 17 sept,2022
//Aim:Coulomb Potential for Hydrogen Atom and Hydrogenic Systems.

clc
clear

for n=1:5
    rmax=n*(1e-9)
    disp("Value of rmax:",rmax)
N=1001
rmin=1*(1e-15)
r=linspace(rmin,rmax,N)
del_r=(r(2)-r(1))
r_new=(r(2:N-1))'

a=0.529D-10            //bohr's radius
Z=1                //atomic no. of hydrogen atom.
ec=9e9                //electrostatic constant = 1/(4*%pi*e0)
h=6.63e-34
m=9.109e-31
e=1.6e-19
C=-((h/(2*%pi))^2)/(2*m*(del_r^2)*e)  //value of the constant
h_cut=(h/(2*%pi))

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
H=KE+V         //Hamilitonian function
[eigfn,eigval]=spec(H)

counter=0
for i=1:N-2
    if eigval(i,i)<0
        counter=counter+1
    end
end
disp("COUNTER",counter)
end

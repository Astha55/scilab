//Date: 16 sept,2022
//Aim: Coulomb Potential for Hydrogen Atom and Hydrogenic Systems

clc
clear

N=1001
rmin=1*(1e-15)
rmax=1*(1e-9) 

r=linspace(rmin,rmax,N)
del_r=(r(2)-r(1))
r_new=(r(2:N-1))'

a=0.529D-10                   //bohr's radius
Z=1
ec=9e9                     //electrostatic constant = 1/(4*%pi*e0)
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
V1=V0(2:N-1)
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

for i=1:counter
    nv=eigfn(:,i)
    norml(i)=(inttrap(r_new,(nv.^2)))      //Normalization of the bound state
    WF(:,i) = eigfn(:,i)./sqrt(norml(i))     //Normalised wavefunction
    exp_r(i)=(inttrap(r_new,(r_new.*(WF(:,i).^2))))    //expectation of r = <r>
    exp_r2(i)=(inttrap(r_new,((r_new.^2).*(WF(:,i).^2))))    //expectation of r^2
    exp_V(i)=(inttrap(r_new,(V1.*(WF(:,i).^2)')))
end

S2 =-(sqrt(-1)*h_cut)
for i=1:N-3
    for j=1:counter
        wf_mid(i,j)= (WF(i,j)+WF(i+1,j))/2
        diff_wf(i,j)= (WF(i+1,j)- WF(i,j))/del_r
        r_mid(i)=(r_new(i)+r_new(i+1))/2
    end
end

for i=1:N-4
    for j=1:counter
        r_mid2(i)= (r_mid(i)+ r_mid(i+1))/2
        wf_mid2(i,j)=(wf_mid(i,j)+wf_mid(i+1,j))/2
        diff2_wf(i,j)= (diff_wf(i+1,j) - diff_wf(i,j))/del_r
    end
end

for i = 1:counter
exp_p(i)=S2*(inttrap(r_mid,(wf_mid(:,i).*diff_wf(:,i))))
exp_p2(i)=-h_cut*h_cut*(inttrap(r_mid2,(wf_mid2(:,i).*diff2_wf(:,i))))
end

exp_KE = (exp_p2/(2*m*e))
var_p =(exp_p2-(exp_p.^2))
sigma_p =sqrt(var_p)
var_r = (exp_r2)-(exp_r.^2)
sigma_r = sqrt(var_r)

u=(sigma_r.*sigma_p)/(h_cut/2)       //Uncertainity Product
TE= exp_KE+exp_V
//disp("Total energy:",TE)

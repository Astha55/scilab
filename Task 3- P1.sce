//Date: 4 Oct,2022
//Aim: Program 1: Eigenvalues and eigenvectors.

clc
clear
clf

x1=1D-15
xn=1D-9
n=2001
x=linspace(x1,xn,n)
del_x=(xn-x1)/(n-1)
disp("The step size is = ",del_x)

e=1.6*10^-19
h=1.054*10^-34
m=9.109*10^-31
C=(1/e)*((-h^2)/(2*m*(del_x^2)))
disp("The C is =",C)

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
xn=x(5:n)
k=9D9
Z=1
s(1)=%inf              //coul
s(2)=0.3D-9            //scr.coul
s(3)=0.5D-9
s(4)=0.7D-9
c(1)='-*b'
c(2)='-^r'
c(3)='-og'
c(4)='-dm'
show_window(1)
for j=1:4
        v(j,:)=-(((k*(e))*(Z.*exp(-x./s(j))))./x)
        vn(j,:)=v(5:n)
        plot(xn,vn(j,:),c(j))
        title("Plot - Coulombic potential and screened coulomb potential ",'fontsize',5,'color','M')
        xlabel('------ x axis -------->','Fontsize',4,'color','black')
        ylabel('------ y axis -------->','Fontsize',4,'color','black')
end
for j=1:4
        v(j,:)=-(((k*(e))*(Z.*exp(-x./s(j))))./x)
        V= diag(v(j,(2:n-1)))
        H=K+V
        [Evec,Eval]=spec(H)
        counter=0
        for i=1:n-2
            if Eval(i,i)<0
                counter=counter+1
            end
        end
        disp("The number of bound states is = ",counter)
        
        disp("The negative Eigen values are = ")
        for i=1:2
            ng(i,1)=Eval(i,i)
        end
        disp(ng)
        ao=0.529*1D-10
        Da=s(j)/ao        // Value of ZD/a0
        disp("The value of ZD/ao is =",Da)
        
        xnew=x(2:n-1)
        xnewt=xnew'
        // Normalization
        for i=1:2
            c=(Evec(:,i).*Evec(:,i))
            in(i)=inttrap(xnew,c)
            NEf(:,i)=(1/sqrt(in(i)))*Evec(:,i)
        end
        show_window(3)
        for i=1:2
            subplot(1,2,i)
            plot(xnewt,NEf(:,i),c(j))
            title("Plot - Normalized Eigenfunctions Vs x for the bound state - "+string(i),'fontsize',5,'color','M')
            xlabel('------ x axis -------->','Fontsize',4,'color','black')
            ylabel('------ y axis -------->','Fontsize',4,'color','black')
            g=gca()
            g.x_ticks=tlist(["ticks","locations","labels"],[5.010D-13;5.000D-10;9.995D-10],["5.010D-13";"5.000D-10";"9.995D-10"]);
            xgrid
        end
end

for j=1:4
        v(j,:)=-(((k*(e))*(Z.*exp(-x./s(j))))./x)
        V= diag(v(j,(2:n-1)))
        H=K+V
        [Evec,Eval]=spec(H)
        counter=0
        for i=1:n-2
            if Eval(i,i)<0
                counter=counter+1
            end
        end
        for i=1:2
            ng(i)=Eval(i,i)
        end
        xnew=x(2:n-1)
        xnewt=xnew'
        // Normalization
        for i=1:2
            c=(Evec(:,i).*Evec(:,i))
            in(i)=inttrap(xnew,c)
            NEf(:,i)=(1/sqrt(in(i)))*Evec(:,i)
            R(:,i)=((NEf(:,i))./xnewt)
        end
        
        for i=1:2
            show_window(4)
            subplot(1,2,i)
            plot(xnewt,R(:,i))
            title("Plot - R(r) Vs x for the bound state - "+string(i),'fontsize',5,'color','r')
            xlabel('------ x axis ------>','Fontsize',4,'color','black')
            ylabel('------ y axis ------>','Fontsize',4,'color','black')
            g=gca()
            g.x_ticks=tlist(["ticks","locations","labels"],[5.010D-13;5.000D-10;9.995D-10],["5.010D-13";"5.000D-10";"9.995D-10"]);
            xgrid
        end
end

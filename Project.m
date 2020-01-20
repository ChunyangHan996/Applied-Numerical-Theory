clc; clear all;
global N D hr A1 A2 A3 b    %Global variables to make the functions work

R=1/100*[6 12 18];       %Given radii in meters
D=4*10^(-7);             %Thermal diffusity
N=60;                    %Number of intervals in r-direction
u0=980*ones(N,1);        %Temperature of the ball when the cooling starts
e=ones(N,1);             %Help-vector for the matrix A
b=zeros(N,1);          %Preallocating the solution vector (RHS)
A2=zeros(N-1,N-1);       %Preallocating the second part of the equation matrix A (RHS)
A2(1,2)=1;               %The second element of the A matrix (A1+A2) must be 2, thus we need to add 1 to A2

k=2;                     %Change k to 1 for radius 6 cm, 2 for radius 12 cm and 3 for radius 18 cm
Tint=[0 100000];                              %A large time interval to get the whole solution (which is then narrowed down)
Ra=R(k);                                      %Selecting the k:th radius
hr=Ra/(N);                                    %Step length in R-direction
A1=spdiags((D/(hr^2))*[e,-2*e,e],-1:1,N-1,N-1); %The "main" equation matrix, first part of A

for i=2:N-2                                    %Forming the other part of the A-matrix
        A2(i,i-1)=-1/i;
        A2(i,i+1)=1/i;
end
A2(N-1,N-2)=-1/(N-1);
A2=(D/(hr^2))*A2;
A31 = [zeros(N-1,1),A1+A2] ;
A32 = (D/(hr^2))*[-6 6 zeros(1,N-2)];
A3 = [A32; A31];
    
for i=1:20                               %Computing 20 times to converge to the lowest Tend
        T=(Tint(1)+Tint(2))/2;
        Tend=T;
        [Tsol,Usol]=ode23s(@(t,u)rhs(t,u,T),[0 Tend], u0(1:N,:)); %Solving for every temperature except at R
        uR=zeros(length(Usol(:,1)),1);   %Preallocating the temperature vector at R
        
        for j=1:length(uR)
            uR(j)=f(Tsol(j),T);          %Computing the temperature at R with the help of the BC  
        end
        
        Usol=[Usol uR];          %Adding the temperature at R 
        dif=abs(diff(Usol,1,2)); %Compute the difference between every column and do so for every row
        dudr = dif/(Ra/N);       %Computing the gradient of u with respect to r
        maxdiff=max(max(dudr,[],2));     %Taking out the max difference in each row then take the max of that 
        
        if maxdiff<6000    %Check if the max value of the gradient is less than 6000. Narrow down the time interval
            Tint(2)=T;      
        else
            Tint(1)=T;
        end
end
e = eig(A3);                %Finding the eigenvalues of A (needs to be <0)

figure(k);
timeh=T/3600           %Time in hours
surf(0:timeh/(length(Usol(:,1))-1):timeh,100*(0:hr:Ra),Usol','EdgeColor','none','LineStyle','none','FaceLighting','phong');
ylabel('Radius, r (cm)');
xlabel('Time, t (hours)');
title(['Temperature u(r,t) for the ball with radius ' num2str(100*Ra) ' cm ']);
       



function rhs=f(t,u,T)

global N D hr A1 A2 A3 b

if t<=T
    temp=980*exp(-3.9*t/T);
else
    temp=980*exp(-3.9);
end

b(N,1)=(D/(hr^2))*(1+(1/(N-1)))*temp;
rhs=(A3)*u+b;

end
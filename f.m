function temp=f(t,T)
if t<T
    temp=980*exp(-3.9*t/T);
else
    temp=980*exp(-3.9);
end
end
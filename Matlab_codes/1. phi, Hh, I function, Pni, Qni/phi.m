function phi= phi(x)
%a= (1+erf(x/sqrt(2)))/2

if x<10 && x>-10, phi=(1+erf(x/sqrt(2)))/2
else phi=0
end
end

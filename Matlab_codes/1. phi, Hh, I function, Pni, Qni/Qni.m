function Qni = Qni(n, i, p, eta1, eta2) 
for j=i:(n-1)
Qni = Sum( Binomial(n, j)*( (1-p)^j )* ( p^(n-j) )*Binomial(n-i-1, j-i)*...
       ( (eta2 /(eta1 + eta2) )^(j-i) ) * ( (eta1 /(eta1 + eta2))^(n-j)))
   
   Qni=(1-p)^n
end
end
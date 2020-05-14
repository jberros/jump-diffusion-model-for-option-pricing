function Pni=Pni(n, i, p, eta1, eta2)
for j=i:(n-1)
Pni= Sum( Binomial(n, j)*( p^j )* ((1-p)^(n-j) )*Binomial(n-i-1, j-i)*...
       ( (eta1/(eta1 +eta2))^(j-i) ) * ( (eta2/(eta1 + eta2))^(n-j)))
      
Pni = p^n
end
end
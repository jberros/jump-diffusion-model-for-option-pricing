function II =II(jj, ll, aa, bb, dd)


if(bb >0 && aa ~=0),
II= -(exp(aa*ll)/aa) * (Table((bb/aa)^(jj-i), {i, 0, jj}) . Table(Hh(i, bb*ll-dd), {i, 0, jj})) +...
((bb /aa)^(jj+1) )*(Sqrt(2*pi)/bb)*Exp(aa*dd/bb + (1/2) * (aa/bb)^2 ) *...
phi(-bb*ll + dd + aa/bb),

elseif(bb<0 && aa <0),
- (Exp(aa*ll)/aa) *(Table((bb/aa)^(jj-i), {i, 0, jj}) . Table(Hh(i, bb*ll-dd), {i, 0, jj}))-...
((bb /aa)^(jj+1) )*(Sqrt(2*Pi)/bb)*Exp(aa*dd/bb + (1/2) * (aa/bb)^2 ) *...
phi(bb*ll - dd - aa/bb),
(bb >0 && aa ==0),
Hh(n+1, bb*ll -dd)/bb
end
end
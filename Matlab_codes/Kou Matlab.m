%% The code to compute the option pricing formulae as in Kou
% A jump diffusion model for option pricing

% Step1. Define phi, Hh, I function, Pni, Qni

function phi= phi(x)
%a= (1+erf(x/sqrt(2)))/2

if x<10 && x>-10, phi=(1+erf(x/sqrt(2)))/2
else phi=0
end
end

function Hh= Hh(n,x)
function temp = temp (n,x)
Hh=sqrt((pi)/2^n)*exp(-(x^2)/2)*(hypergeom((n+1)/2,1/2,(x^2)/2)/(sqrt(2)*gamma(1+n/2)))
-x*(hypergeom((n/2)+1/2,3/2,(x^2)/2)/(gamma((1+n)/2)))

if x>=-6 && x<10, Hh=(1/factorial(n))*integrate((t-x)^n*exp(-t^2/2),t,x,infinity)
else Hh=0,
    temp =(x+sqrt(x^2+4*n)/2);
    Hh=(1/factorial(n))*(integrate((t-x)^n*exp(-t^2/2),t,x,temp-3)+...
integrate((t-x)^n*exp(-t^2/2),t,temp-3,temp-1)+...
integrate((t-x)^n*exp(-t^2/2),t,temp-1,temp)+...
integrate((t-x)^n*exp(-t^2/2),t,temp,temp+1))+...
integrate((t-x)^n*exp(-t^2/2),t,temp+1,temp+3)+...
integrate((t-x)^n*exp(-t^2/2),t,temp+3,infinity)
end
end
end

function Ha = Ha(m,y)
Ha = (-y/n)*Ha(n - 1) + (1/n)*Ha(n - 2); 
Ha(-1) = exp(-y^2/2); 
Ha(0) = sqrt(2*pi)*phi(-y);
Ha(m)
end

function II =II(jj, ll, aa, bb, dd)


if(bb >0 && aa ~=0),
II=(exp(aa*ll)/aa) * (Table((bb/aa)^(jj-i), {i, 0, jj}) . Table(Hh(i, bb*ll-dd), {i, 0, jj})) +...
((bb /aa)^(jj+1) )*(Sqrt(2*pi)/bb)*Exp(aa*dd/bb + (1/2) * (aa/bb)^2 ) *...
phi(-bb*ll + dd + aa/bb),

elseif(bb<0 && aa <0),
- (Exp(aa*ll)/aa) *(Table((bb/aa)^(jj-i), {i, 0, jj}) . Table(Hh(i, bb*ll-dd), {i, 0, jj})) -...
((bb /aa)^(jj+1) )*(Sqrt(2*Pi)/bb)*Exp(aa*dd/bb + (1/2) * (aa/bb)^2 ) *...
phi(bb*ll - dd - aa/bb),
(bb >0 && aa ==0),
Hh(n+1, bb*ll -dd)/bb
end
end

function Pni=Pni(n, i, p, eta1, eta2)
for j=i:(n-1)
Pni= Sum( Binomial(n, j)*( p^j )* ((1-p)^(n-j) )*Binomial(n-i-1, j-i)*...
       ( (eta1/(eta1 +eta2))^(j-i) ) * ( (eta2/(eta1 + eta2))^(n-j)))
      
Pni = p^n
end
end

function Qni = Qni(n, i, p, eta1, eta2) 
for j=i:(n-1)
Qni = Sum( Binomial(n, j)*( (1-p)^j )* ( p^(n-j) )*Binomial(n-i-1, j-i)*...
       ( (eta2 /(eta1 + eta2) )^(j-i) ) * ( (eta1 /(eta1 + eta2))^(n-j)))
   
   Qni=(1-p)^n
end
end

%Step 2a. calculate the call price 
 
function cprob= cprob(mu, eta1, eta2, la, p, sig, aa,  bigT, nStep)
    function IITwo = II(jj, ll, aa, bb, dd)
        for k=1:nStep
    IITwo = Table(  II(k-1, aa - mu * bigT, -eta1, -1/(sig*Sqrt(bigT)), -(sig*Sqrt(bigT))*eta1));
        end
    end

    function IIFour = II(jj, ll, aa, bb, dd)
        for k=1:nStep
    IIFour = Table( II(k-1, aa - mu * bigT, eta2, 1/(sig*Sqrt(bigT)), -(sig*Sqrt(bigT))*eta2));
        end
    end

    function PiN=PiN(n) 
        PiN= Exp(-la*bigT)*((la*bigT)^n) /(factorial(n));
    end

    function PiNPni
        for n=1:nStep
            for k=1:n
  PiNPni = Table (PiN(n) * Pni(n, k, p, eta1, eta2) * ((sig*Sqrt(bigT)*eta1)^k));
            end
        end
    end

    function PiNQni
        for n=1:nStep
            for k=1:n
  PiNQni = Table (PiN(n) * Qni(n, k, p, eta1, eta2) * ((sig*Sqrt(bigT)*eta2)^k));
            end
        end
    end

                 
    function sec
        for n=1:nStep
            for k=1:n
   sec = Sum(PiNPni((n,k)) * IITwo((k)));
            end
        end
    end

    function fourth
        for n=1:nStep
            for k=1:n
  fourth = Sum(PiNQni ((n, k)) * IIFour((k)));
            end
        end
    end

cprob=(sec * Exp(((sig*eta1)^2)*bigT/2) + fourth * Exp(((sig*eta2)^2)*bigT/2)) /( Sqrt(2*Pi) * sig * Sqrt(bigT) ) +... 
    Exp(- la* bigT) * phi (-(aa- mu*bigT)/(sig*Sqrt(bigT)))
end

% The function callOR computes the regular call price for stock options

function callOR = callOR(eta1, eta2, la, p, sig, rr, bigS, bigK,  bigT, nStep)

    function zetaaOR
        zetaaOR = p*eta1 /(eta1 -1) + (1-p)*eta2/(eta2+1) -1;
    end

    function tempaa1OR
      tempaa1OR =  rr + sig*sig/2 - la*zetaaOR;
    end

    function tempaa2OR
      tempaa2OR = tempaa1OR - sig*sig;
    end

callOR = bigS *...
 cprob(tempaa1OR, eta1 - 1, eta2 + 1, la *(1+zetaaOR), p*eta1/((1+zetaaOR)*(eta1-1)), sig, Log(bigK/bigS),  bigT, nStep) -...
  bigK * Exp(- rr*bigT) *cprob(tempaa2OR, eta1, eta2, la, p, sig, Log(bigK/bigS),  bigT, nStep)
end

% The function call computes the regular call price for futures options

function call = call(eta1, eta2, la, p, sig, bond, bigF, bigK,  bigT, nStep)

    function zetaa
      zetaa = p*eta1 /(eta1 -1) + (1-p)*eta2/(eta2+1) -1;
    end

    function tempaa1
      tempaa1 =   sig*sig/2 - la*zetaa;
    end

    function tempaa2
      tempaa2 = tempaa1 - sig*sig;
    end

   call = bond* ( bigF *...
 cprob(tempaa1, eta1 - 1, eta2 + 1, la *(1+zetaa), p*eta1/((1+zetaa)*(eta1-1)), sig, Log(bigK/(bigF)),  bigT, nStep) -...
  bigK *...
  cprob(tempaa2, eta1, eta2, la, p, sig, Log(bigK/(bigF)),  bigT, nStep))
end
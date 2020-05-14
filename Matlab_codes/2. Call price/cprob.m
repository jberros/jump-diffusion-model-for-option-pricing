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
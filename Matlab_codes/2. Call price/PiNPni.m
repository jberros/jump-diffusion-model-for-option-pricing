function PiNPni
        for n=1:nStep
            for k=1:n
  PiNPni = Table (PiN(n) * Pni(n, k, p, eta1, eta2) * ((sig*Sqrt(bigT)*eta1)^k));
            end
        end
    end

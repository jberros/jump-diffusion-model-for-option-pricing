function PiNQni
        for n=1:nStep
            for k=1:n
  PiNQni = Table (PiN(n) * Qni(n, k, p, eta1, eta2) * ((sig*Sqrt(bigT)*eta2)^k));
            end
        end
    end
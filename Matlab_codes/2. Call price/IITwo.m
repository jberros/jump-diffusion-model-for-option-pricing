function IITwo = IITwo(jj, ll, aa, bb, dd)
        for k=1:nStep
    IITwo = Table(  II(k-1, aa - mu * bigT, -eta1, -1/(sig*Sqrt(bigT)), -(sig*Sqrt(bigT))*eta1));
        end
    end
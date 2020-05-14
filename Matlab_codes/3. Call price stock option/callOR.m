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
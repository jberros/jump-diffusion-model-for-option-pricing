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
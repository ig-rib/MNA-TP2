function nextX = symmetricIntegrator(q, h, X, gammas, k)
  n = floor(q/2);
  spmd(q)
     for j = 1:n
        if labindex==j
            for s = 1:labindex
                X = lieTrotterPlus(h/labindex, X, k);
            end
        end
     X = gammas(labindex)*X;
         if labindex == j+n
             for s = 1:labindex-n
                X = lieTrotterMinus(h/labindex, X, k);
             end
         X = gammas(labindex-n)*X;    
         end
     end
  end
  for s = 2:q
    x{1} = x{1} + x{s};
  end
  nextX = x{1};
end
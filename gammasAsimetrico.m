function gammas = gammasAsimetrico(q) 
    M=ones(q);
    
    for i=2:q
        for j=1:q
            M(i,j) = j^(-i+1);
        end
    end
    
    gammas = M\([1 zeros(1, q-1)]');
end
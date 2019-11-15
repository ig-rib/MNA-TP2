function gammas = gammasSimetrico(q) 
    

    n=floor(q/2);    
    M=ones(n);
    
    for i=2:n
        for j=1:n
            M(i,j) = j^(2*(-i+1));
        end
    end
    
    gammas = M\([1/2 zeros(1, n-1)]');
end
function Z = Z_indep_Ising( exph_l )
%Z_INDEP_ISING 
% compute normalizing coefficient for Ising model with independent neurons

Z = [1 exph_l(1)];

for n = 2:length(exph_l)
    Z = conv(Z,[1 exph_l(n)]);
end

Z = sum(Z);

end


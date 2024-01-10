function [mantissa, exponent] = get_scientific_decomposition(x)
    mantissa = x(:).*10.^ceil(-log10(abs(x(:)+(x==0))));
    exponent = int8(floor(log10(abs(x(:)+(x==0)))));
end
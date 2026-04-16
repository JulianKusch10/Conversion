function [Observ_julia, Observ_mat] = get_arrays(Observ_julia, Observ_mat)
    n_julia = length(Observ_julia); 
    Observ_julia = Observ_julia(2:end); 
    Observ_mat = Observ_mat(1:n_julia-1); 
    Observ_mat = Observ_mat(:); 
    Observ_julia = Observ_julia(:);
end
function HTherm(Params, psi)
    T = Params.T 
    as = Params.as/Params.a0 

    coeffs_A = h5read("HThermFit_data.h5", "/coeffs_A")
    coeffs_B = h5read("HThermFit_data.h5", "/coeffs_B")

    n = abs.(psi).^2
    # n = 1000

    function A(as, T)
        coeffs_A[1] + coeffs_A[2]*as + coeffs_A[3]*T + coeffs_A[4]*as^2 + coeffs_A[5]*as*T + coeffs_A[6]*T^2
    end
    function B(as, T)
        coeffs_B[1]/(coeffs_B[2]*T + coeffs_B[3]*sqrt(T) + coeffs_B[4]) + coeffs_B[5]*exp(-coeffs_B[6]*as)
    end
    function HT_model(as, T, n)
        A(as, T)*exp.(-B(as, T)*sqrt.(n))
    end
    HT = HT_model(as, T, n)
    return HT
end
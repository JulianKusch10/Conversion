using Interpolations
using NumericalIntegration
using Dierckx

function EthermInterp(Params, psi)

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

    logdenslist = range(-9, stop=6, length=1000)   # 1e3 → 1000 points
    densList = 10 .^ logdenslist  
    densList = [0; densList]
    
    HT = HT_model(as, T, densList)
    
    EthTable = cumul_integrate(densList, HT)
    # itp = interpolate((densList,), EthTable,Gridded(Cubic(Line(OnGrid()))))
    # Ftherm = extrapolate(itp, Line())   # linear extrapolation outside range
    Ftherm = Spline1D(densList, EthTable, k=3, s=0, bc="extrapolate")  # cubic spline
    return Ftherm
end
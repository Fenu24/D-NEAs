function [name,a,e,inc,H,H_sigma,D,D_sigma,dadt,dadt_sigma,...
    P,P_sigma,Q,SNR, ref] =  tableToVector(T,valSNR,vale)

    index = (T.SNR > valSNR & T.e < vale);
    
    name       = T.Name(index);
    a          = T.a(index);
    e          = T.e(index);
    inc        = T.i(index);
    H          = T.H(index);
    H_sigma    = T.H_std(index);
    D          = T.D(index)*1000;
    D_sigma    = T.D_std(index)*1000;
    dadt       = T.dadt(index);
    dadt_sigma = T.dadt_std(index);
    P          = T.P(index);
    P_sigma    = T.P_std(index);
    Q          = T.Q(index);
    SNR        = T.SNR(index);
    ref        = T.Ref(index);
end
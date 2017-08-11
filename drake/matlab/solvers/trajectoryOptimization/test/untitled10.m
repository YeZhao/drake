
%%
switch pp2
    case 1 % x0
        dL = sum(dfdSig2*dfdSig1*dfdx0,1);
        dL = dL(:);
    case 2 % u0
        dL = sum(dfdSig2*dfdSig1*dfdu0,1);
        dL = dL(:);
    case 3 % sigma0
        dL = sum(dfdSig2*dfdSig1*dfdSig0,1);
        dL = dL(:);
    case 4 % x1
        dL = sum(dfdSig2*dfdx1,1);
    case 5 % u1
        dL = sum(dfdSig2*dfdu1,1);
    case 6 % x2
        dL = sum(dfdx2,1)';
    case 7 % u2
        dL = sum(dfdu2(:));
end

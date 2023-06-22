function SOP = SOP_NakagamiM_N(beta, eta, OmegaB, OmegaE, mB, mE ,choice)
    %t = 2*sigma^2;
    switch choice
    case 1
        % SOP = int_0^{infty} F_{gAB}(beta + eta*x)f_{gAE}(x)dx
        
        term1 =  (exp(-beta./OmegaB)) ./ ( gamma(mE).*( (1+ (OmegaE.*eta./OmegaB) ).^mE ) ) ;
        sumT = 0;
        for n = 0:mB-1
    %         term2 = ((beta./OmegaB).^n)*(1/(factorial(n)));
            term3 = 0;
            for k = 0:n
                term3 = term3 + (  nchoosek(n,k)*(( eta./( (beta./OmegaE).*(1+ (OmegaE.*eta./OmegaB) ) ) ).^k)*gamma(mE+k)  );

                %term1 = nchoosek(n,ln) * ( gamma(k+ln) ./ ( factorial(n)*gamma(k)*(tAB.^n)*(tAE.^k) ) ) .* (2.^(Rs.*ln)).*( (2.^Rs - 1).^(n-ln) );
                %term2 = ( (bJ.^ln)./(aJ.^n) ) .* ( ( (2.^Rs)*(bJ./aJ).*(1./tAB) + (1./tAE) ).^(-k-ln) ) .* exp( -(1./tAB).*(2.^Rs - 1).*(1./aJ) );
                %sumT =  sumT + term1.*term2;
            end
    %         sumT = sumT + term2.*term3;
            sumT = sumT + ((beta./OmegaB).^n).*(term3/(gamma(n+1)));

        end
        SOP = 1 - term1.*sumT;
        % Approximation: SOP = F_{gAB}(beta)
        %SOP = (1./(gamma(mB)))*gammainc(mB,beta./OmegaB);
    end
    SOP(SOP<0)=0;
end

function [mu0, mu1, etaVect] = CPWL_param(aVect, bVect, a0)
    
    % N=length(aVect);
    % mVectTilde=[m1; m2;.. mj;...mN-1]; with j=1,...,N-1 
    mVectTilde = (bVect(2:end)-bVect(1:end-1))./(aVect(2:end)-aVect(1:end-1));    
    m0 = mVectTilde(1);
    mN = mVectTilde(end);
    mVect = [m0 mVectTilde mN];
      
    %nuVect = 0.5*(bVect(2:end) - bVect(1:end-1));
    etaVect = 0.5*(mVect(2:end)-mVect(1:end-1));
    
  
    mu1 = 0.5*(m0+mN);    
    mu0 = a0 - etaVect*abs(aVect'); %- nuVect'* sign(aVect);

end




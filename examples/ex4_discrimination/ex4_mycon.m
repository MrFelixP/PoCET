function [c,ceq] = ex4_mycon(u,pcesysH,pcesysM,simoptions)
 simH = PoCETsimGalerkin(pcesysH,'ex4_ODE_H',[],simoptions,'u_v',u); % simulate system H
 momH = PoCETcalcMoments(pcesysH,pcesysH.MomMats,simH.x_2.pcvals(:,end)); % calc. moments
 betaH = calcBeta4(momH); % fit 4-parameter beta distribution

 simM = PoCETsimGalerkin(pcesysM,'ex4_ODE_M',[],simoptions,'u_v',u); % simulate system M
 momM = PoCETcalcMoments(pcesysM,pcesysM.MomMats,simM.x_2.pcvals(:,end)); % calc. moments
 betaM = calcBeta4(momM); % fit 4-parameter beta distribution

 c = calcMDCbeta(betaH,betaM); % use Hellinger distance as measure of PDF overlap
 ceq = []; % no equality constraints


 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 %% 
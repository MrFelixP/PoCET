function PSImap = get_PSImap(N_xi,P)
% PoCET internal function. Do not call directly!

PSI = zeros(1,N_xi);
P_tilde = factorial(N_xi+P)/(factorial(N_xi)*factorial(P));
PSImap = zeros(P_tilde,N_xi);
for i = 2:P_tilde
 PSI = next_PSI_k(PSI,P);
 PSImap(i,:) = PSI;
end
end

function PSI = next_PSI_k(PSI,Pmax)
Pcur = sum(PSI);                 % get current polynomial order
if Pcur <= Pmax && PSI(end)<Pmax % if current order < maximum order ...
 if numel(PSI) > 1
  if PSI(end)==Pcur               % ... and xi_N is at current order ...
   PSI(1) = Pcur+1;               % ... go to next higher order ...
   PSI(end) = 0;                  % ... and reset xi_N
  else
   n0 = find(PSI,1);              % find first non-zero element
   if n0 ~= 1                     % if first non-zero element is not equal to
    PSI(1) = PSI(n0) - 1;
    PSI(n0+1) = PSI(n0+1) + 1;
    PSI(n0) = 0;
   else
    PSI(1) = PSI(1) - 1;
    PSI(2) = PSI(2) + 1;
   end
  end
 else
  PSI = PSI + 1;
 end
end
end
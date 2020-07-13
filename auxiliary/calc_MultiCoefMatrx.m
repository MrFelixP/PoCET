function [C,MATS] = calc_MultiCoefMatrx(PoCETsystem,xis,pars,order,PHI_rj,wj)
% PoCET internal function. Do not call directly!

[~,sIDX] = sort(xis);
xis = xis(sIDX);      % sort input
pars = pars(sIDX);    % sort input

opts = PoCETsystem.pce.options;
n_xis = numel(xis); 
if order == 0
 C = zeros(opts.n_phi,1);
 tmpCOEF = 1;
 for i = 1:numel(pars) % Liste der Mischterme abarbeiten
  tmpCOEF = tmpCOEF * PoCETsystem.parameters(pars(i)).pce(1);  % zugehoerige Koeffizienten berechnen sys.p(pars(j)).co(coefnums(i,j)) 
  if xis(i) ~= 0
   C(xis(i)+1) = PoCETsystem.parameters(pars(i)).pce(xis(i)+1);
  end
 end
 C(1) = tmpCOEF; % in Gesamtmatrix speichern
 MATS = {[];[];C};
 
else
 if n_xis == 0 || all(xis == 0)
  C = calc_CoefMatrx(opts,1,order,PHI_rj,wj);
  tmpCOEF = 1;
  for i = 1:numel(pars) % Liste der Mischterme abarbeiten
   tmpCOEF = tmpCOEF * PoCETsystem.parameters(pars(i)).pce(1);  % zugehoerige Koeffizienten berechnen sys.p(pars(j)).co(coefnums(i,j)) 
  end
  C = C * tmpCOEF; % in Gesamtmatrix speichern
  MATS = {[];[];C};
 else
  % store all nonzero pce coefficients for each parameter
  phis = cell(1,n_xis); 
  for i = 1:n_xis
   phis{i} = find(PoCETsystem.parameters(pars(i)).pce)-1;
  end
  
  % create combinations of all parameter-coefficients
  allmats = xicombs(phis{:});
  
  % calculate the matrix
  allmats = allmats + 1;
  MATS  = cell(3,size(allmats,1));
  C = sparse(opts.n_phi,opts.n_phi^order); % preallocation
  for i = 1:size(allmats,1) % cycle through list of combinations
%    keyboard
   tmpCOEF = 1;
   tmpMAT = calc_CoefMatrx(opts,allmats(i,:),order,PHI_rj,wj); % calculate matrix
   MATS{1,i} = pars;
   MATS{2,i} = allmats(i,:);
   MATS{3,i} = tmpMAT;
   for j = 1:n_xis
    tmpCOEF = tmpCOEF * PoCETsystem.parameters(pars(j)).pce(allmats(i,j));  % calculate coefficient
   end
   C = C + tmpCOEF * tmpMAT; % put it all together
  end
 end
end

function A = xicombs(varargin)
 %error(nargchk(1,Inf,nargin)) ;
 narginchk(1,Inf);
 NC = nargin;
 ii = NC:-1:1;
 if NC > 1
  [A{ii}] = ndgrid(varargin{ii}); % flip using ii
  A = reshape(cat(NC+1,A{:}),[],NC); % concatenate
 elseif NC==1,
  A = varargin{1}(:); % nothing to combine
 else % NC==0
  A = zeros(0,0); % nothing
 end
end

end

%% OLD CODE (REMOVED BY 05/30/15)
% Liste 'phis' aller -notwendigen- univariaten Polynome erstellen
% Liste -aller- (multi- und univariaten) Polynome erstellen
%   PSImap = get_PSImap(opts.n_xi,opts.order); 
% add first polynomial -> mean will always be set!
%   phis = 0;
% zudem Liste fuer Polynome 2. Ordnung vorbereiten
%   phi2 = [];
%   keyboard
% Anzahl der PCE-Koeffizienten pro Parameter ermitteln
% -> gibt notwenige Ordnungen der Polynome vor
% -> momentan maximal 3 Koeffizienten (Polynome der Ordnungen 0, 1 & 2)
%   n_coefs = zeros(1,n_xis); 
%   for i = 1:numel(n_coefs) 
%    n_coefs(i) = nnz(PoCETsystem.parameters(pars(i)).pce);
%    % phis = find(PoCETsystem.parameters(pars(i)).pce)-1;
%    if n_coefs(i) == 2     % falls PCE-Koeffizient 1. Ordnung vorhanden...
%     tmp1 = find(PSImap(:,xis(i))==1,1,'first')-1;
%     phis = [phis tmp1]; % entsprechendes Polynom zur Liste hinzufuegen
%    elseif n_coefs(i) == 3 % falls PCE-Koeffizienten 1. & 2. Ordnung vorhanden...
%     tmp2 = find(PSImap(:,xis(i))==2,1,'first')-1; % in Liste aller Polynome nach dem univariaten Polynom 2. Ordnung suchen...
%     phis = [phis xis(i) tmp2];  % entsprechende Polynome zur Liste hinzufuegen...
%     phi2 = [phi2; xis(i) tmp2]; % und merken, dass sie zusammengehoeren
%    end
%   end

% Liste aller Mischterme verschiedener Ordnungen erstellen
%      allmats = unique(sort(nchoosek(phis,n_xis),2),'rows'); % *MAGIC*
%   keyboard
%   if n_xis > 1
%    allmats = [zeros(1,n_xis); allmats];
%   end
% Mischterme entfernen, die Polynome erster und zweiter Ordnung der
% gleichen Zufallsvariable enthalten -> treten beim Ausmultiplizieren
% nicht auf, werden hier aber erstmal mit erstellt (deshalb oben merken,
% welche Polynome zusammengehoeren)
%   for i = 1:size(phi2)
%    for j = 1:size(allmats,1)
%     if ismember(phi2(i,1),allmats(j,:)) && ismember(phi2(i,2),allmats(j,:))
%      allmats(j,:) = 0; % entsprechende Terme erstmal durch Nullzeilen ersetzen
%     end
%    end
%   end
%   allmats = unique(allmats,'rows'); % alle Nullzeilen (ausser die erste)
%   wieder entfernen
% Liste aller zu den Mischtermen gehoerender Koeffizienten in richtiger
% Reihenfolge erstellen (vorgegeben durch die Reihenfolge in 'xis')
%   keyboard
%   coefnums = ones(size(allmats));
%   for i = 1:size(allmats,1) % Liste aller Mischterme durchgehen
%    for j = find(allmats(i,:)); % Polynome mit Ordnungen > 0 suchen
%     xiid = find(xis==allmats(i,j)); % #Polynome 1. Ordnung entsprechen #Xi
%     if ~isempty(xiid)
%      coefnums(i,xiid) = 2; % 2. Koeffizient gehoert zu Polynom 1. Ordnung
%     else % Falls #Polynom nicht i 'xis' enthalten ist, muss es 2. Ordnung sein
%      [col,~] = find(phi2==allmats(i,j)); % entsprechendes Xi suchen
%      coefnums(i,xis==phi2(col,1)) = 3; % 3. Koeffizient gehoert zu Polynom 2. Ordnung
%     end
%    end
%   end
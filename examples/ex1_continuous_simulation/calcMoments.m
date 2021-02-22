function moments = calcMoments(MomMats,xi,type,M)
steps = size(xi,2);
raw = zeros(M,steps);
raw(1,:) = xi(1,:).*MomMats(1).val;
if M >= 2
 raw(2,:) = sum(xi.^2.*(MomMats(2).val*ones(1,steps)),1);
end
if M >= 3
 for iM = 3:M
  val = MomMats(iM).val;
  phi = MomMats(iM).phi;
  for k = 1:steps
   raw(iM,k) = val'*prod(reshape(xi(phi,k),[size(phi),1]),2);
  end
 end
end

moments = zeros(M,steps);
switch lower(type)
 case 'raw'
  moments = raw;
 case 'central'
  moments(1,:) = raw(1,:); 
  if M >= 2; moments(2,:) = raw(2,:)-raw(1,:).^2; end
  if M >= 3; moments(3,:) = ...
     (2*raw(1,:).^3 - 3*raw(1,:).*raw(2,:) + raw(3,:))./sqrt(raw(2,:)-raw(1,:).^2).^3; end
  if M >= 4; moments(4,:) = ...
     (-3*raw(1,:).^4 + 6*raw(1,:).^2.*raw(2,:) - 4*raw(1,:).*raw(3,:) + raw(4,:))./(raw(2,:)-raw(1,:).^2).^2 - 3; end
end
end
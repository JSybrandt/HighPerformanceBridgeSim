function [KTotal] = KInsert(ElementK,DOFLocations,TNDOF)

KTotal = zeros(TNDOF);

for i=1:length(DOFLocations)
   for j=1:length(DOFLocations)
        KTotal(DOFLocations(i),DOFLocations(j)) = ElementK(i,j); 
   end
end

end
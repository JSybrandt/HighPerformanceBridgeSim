function [KTotal] = KInsert2(ElementK,DOFLocations,TNDOF)

KTotal = zeros(TNDOF,1);

for i=1:length(DOFLocations)
        KTotal(DOFLocations(i),1) = ElementK(i,1); 
end

end
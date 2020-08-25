function J = JAC_simple( t, y )

global B1 B2 I1 I2 k_hexanal_transfer k_pentanal_transfer hexanal pentanal
o2 = 2;

J =  B1  +  B2 * ( sparse(1:size(I2,1),1:size(I2,1), I2 * y )* I1 + sparse(1:size(I1,1),1:size(I1,1), I1 * y ) * I2 );

J( o2,o2 ) = 0;
J( hexanal,hexanal ) =  J( hexanal,hexanal ) - k_hexanal_transfer;
J( pentanal,pentanal ) =  J( pentanal,pentanal ) - k_pentanal_transfer;


end

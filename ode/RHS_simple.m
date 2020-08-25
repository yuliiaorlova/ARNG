
function dy = RHS_simple( t, y )
o2 = 2;

global B1 B2 id  k_hexanal_transfer k_pentanal_transfer  c_hexanal_equil c_pentanal_equil hexanal pentanal ind_term alpha diff


dy =  B1 * y +    B2 * prod( y( id ), 2 );
%     pos_flux(count, :) = B1_pos * y + B2_pos * prod( y( id ), 2 );
%     neg_flux(count, :) = B1_neg * y + B2_neg * prod( y( id ), 2 );


dy( o2 ) = 0;
dy(hexanal)  = dy(hexanal) +  k_hexanal_transfer * ( c_hexanal_equil - y(hexanal) );
dy(pentanal)  = dy(pentanal)+ k_pentanal_transfer * ( c_pentanal_equil - y(pentanal) );


end


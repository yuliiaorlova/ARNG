function [O2_upt, PV, p_y, hex_sum, pent_sum, y2, t2] = model_run_1( )
%ODE_SOLVER input arduments: rate coefficients, 
% output - oxygen uptake, peroxide value, lumped concentrations, hexanal,
% pentanal, all concentration profiles and time steps
global A Reactions Molecules weight type y0 hexanal pentanal oxygen_containing_species 
global B1 B2 id g I1 I2 ind_term alpha c_hexanal_equil c_pentanal_equil k_hexanal_transfer k_pentanal_transfer
global B1_pos B1_neg B2_pos B2_neg 
global count 
count = 1;


c_ba_0= 3.7106;

tic
k_r = weight(Reactions);
S  = Molecules;
M1 = sum( A( S, Reactions ) ) == 1;
M2 = sum( A( S, Reactions ) ) == 2;
y0 = weight(Molecules);

R1 = Reactions( M1 );
R2 = Reactions( M2 );
type_2 = type(R2);



K1 = k_r(M1);
K2 = k_r(M2);
K1_diag = spdiags(k_r( M1 ), 0, length(k_r( M1 )), length(k_r( M1 )));
K2_diag = spdiags(k_r( M2 ), 0, length(k_r( M2 )), length(k_r( M2 )));
B1 = sparse( A( R1, S )' * K1_diag * A( S, R1 )'  -  diag( A( S, R1 ) * K1 ) );

% indexes of second order reactants, different species
ASR2  = A( S, R2 );
ASR2s = ASR2;
ASR2s( :, sum( ASR2 > 0 ) == 1  ) = 0;
id = 0;
[ a b ] = find( ASR2s );

% indexes of second order reactants, same species
id( b( 1 : 2 : end ), 1 : 2 ) = [ a( 1 : 2 : end ) a( 2 : 2 : end ) ];
ASR2s = ASR2;
ASR2s( :, sum( ASR2s > 0 ) == 2  ) = 0;
[ a b ] = find( ASR2s );
id( b, : ) = [ a  a ];

B2 = sparse( ( A( R2, S )' - A( S, R2 ) ) * K2_diag );

% if size(B2,2)<length(Molecules)
%     B2(length(Molecules), length(Molecules)) = 0;
% end;

I1 = sparse( 1:length(id(:,1)), id(:,1), 1 );

if size(I1,2)<length(Molecules)
   I1(size(I1,1), length(Molecules)) = 0;

end;

 I2 = sparse( 1:length(id(:,2)), id(:,2), 1 );
 
if size(I2,2)<length(Molecules)
  I2(size(I2,1), length(Molecules)) = 0;
end;

options = odeset( 'reltol', 1e-5, 'abstol', 1e-14, 'Stats', 'off' );
options.OutputFcn = @ShowStatus;
options.Jacobian=@JAC_simple;

B1_pos = B1;
B1_neg = B1;
B2_pos = B2;
B2_neg = B2;

B1_pos(B1_pos<0) = 0;
B2_pos(B2_pos<0) = 0;
B1_neg(B1_pos>0) = 0;
B2_neg(B2_pos>0) = 0;

g.progress = 0;
g.t_end = 1000 * 3600;
warning off;
prep_time = toc;
tic
disp( [ 13 '[           ' ] )

[ t, y1 ] = ode15s( @( t, y ) RHS_simple( t, y ), [ 0 g.t_end ], y0,  options ); 

sol_time = toc;
disp( [ 'Preparation time ' num2str( prep_time ) ] )
disp( [ 'Integration time ' num2str( sol_time ) ] )

t2 = t/3600;
rr = find(t2>.01);
t2 = t2(rr);
y2 = y1(rr, :);



 hex_sum  = cumtrapz(t2(:), -k_hexanal_transfer * ( c_hexanal_equil - y2(:, hexanal) ));
 pent_sum  = cumtrapz(t2(:), -k_pentanal_transfer * ( c_pentanal_equil - y2(:, pentanal)));
 
 hex_sum  = hex_sum/c_ba_0;
 pent_sum  = pent_sum/c_ba_0;


% on each timestep multiply concentrations with oxygen_containing_species
ox_species_concentration = zeros(length(t2), 1);

for i = 1 : length(t2)
    ox_species_concentration(i) = sum(y2(i,:)'.*oxygen_containing_species) - y2(i,2);
end

O2_upt = ox_species_concentration/c_ba_0;


peroxides = [2, 3, 9, 10, 14, 21, 22, 25, 27, 32, 39];
crosslinks = [14, 15, 16, 26, 27, 28, 39];

[ p_y ] = lumping_procedure( y2 );
p_y = p_y/c_ba_0;
p_y(:, crosslinks) = p_y(:, crosslinks)/2;
PV = sum(p_y(:, peroxides), 2);


end

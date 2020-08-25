function [  ] = rules(  )
global pattern;
global pattern_non_reactive;
%structures for different types of reactions
global R11 R12 R21 R22 R23
global ks 
% first order rules (1 -> 1)

% 11.1 initiation
R11(1).reactant1 = pattern(3); % Co(3)OOH
R11(1).product1 = pattern(1); % Co(2)
R11(1).k = 0;
R11(1).type = 2; 


% first order rules (1 -> 2)
% 12.1 cobaltperoxide decomposition

R12(1).reactant1 = pattern(25); % ROOCo
R12(1).product1 = pattern_non_reactive(23); % Co(II) + ROO.
R12(1).k = 0;
R12(1).type = 4; 

% 12.2 cobaltperoxide decomposition
R12(end+1).reactant1 = pattern(32); % R(=O)OOCo
R12(end).product1 = pattern_non_reactive(24); % Co(II) + R(=O)OO.
R12(end).k = 0;
R12(end).type = 4;

% 12.3 hydroperoxide decomposition
R12(end+1).reactant1 = pattern(10); % ROOH
R12(end).product1 = pattern_non_reactive(3); % RO. + .OH
R12(end).k = 0;
R12(end).type = 5; 

% 12.4 hydroperoxide decomposition
R12(end+1).reactant1 = pattern(31); % ROOH
R12(end).product1 = pattern_non_reactive(18); % RO. + .OH
R12(end).k = 0;
R12(end).type = 5; 

% 12.5 beta scission: conjugated alkoxy decomposes into dienal and alkyl
% radical D1
R12(end+1).reactant1 = pattern(17); % conjugated alkoxy radical
R12(end).product1 = pattern_non_reactive(4); % dienal and alkyl radical
R12(end).k = 0;
R12(end).type = 6; 


% second order rules (2 -> 1)

% 21.1 Oxidation of catalyst
R21(1).reactant1 = pattern(1); % Co(II)
R21(1).reactant2 = pattern(4); % O2 
R21(1).product1 = pattern(2); % Co(III)OO.
R21(1).k = 0;
R21(1).type = 9;

% 21.2 Oxidation 
R21(end+1).reactant1 = pattern(8); % R.
R21(end).reactant2 = pattern(4); % O2 
R21(end).product1 = pattern_non_reactive(10); % ROO.
R21(end).k = 0;
R21(end).type = 10; 

% 21.3 Oxidation of radical after beta-scission 
R21(end+1).reactant1 = pattern(20); % RC(O).
R21(end).reactant2 = pattern(4); % O2
R21(end).product1 = pattern(21); % RC(O)OO.
R21(end).k = 0;
R21(end).type = 10; 

% 21.4 oxidation of radical after beta-scission 
R21(end+1).reactant1 = pattern(19); % R.
R21(end).reactant2 = pattern(4); % O2
R21(end).product1 = pattern_non_reactive(7); % ROO.
R21(end).k = 0;
R21(end).type = 10;

% 21.5 oxidation of radical after Beta scission D2
R21(end+1).reactant1 = pattern(39); % 
R21(end).reactant2 = pattern(4); % O2
R21(end).product1 = pattern_non_reactive(26); % R_OO._cl_cl
R21(end).k = 0;
R21(end).type = 10;

% second order rules (2 -> 2)
%__________________________________

% 22.1 Initiation with Co 
R22(1).reactant1 = pattern(2); % Co(III)OO.
R22(1).reactant2 = pattern(7); % RH
R22(1).product1 = pattern(3); % Co(III)OOH
R22(1).product2 = pattern_non_reactive(1); % R.
R22(1).k = 0;
R22(1).type = 11;


% 22.2 Initiation with Co 
R22(end+1).reactant1 = pattern(2); % Co(III)OO.
R22(end).reactant2 = pattern(7); % RH
R22(end).product1 = pattern(3); % Co(III)OOH
R22(end).product2 = pattern_non_reactive(2); % R.
R22(end).k = 0;
R22(end).type = 11;


% 22.3 Hydrogen abstraction from C_BA by OO. k_H_dc
R22(end+1).reactant1 = pattern(9); % ROO.
R22(end).reactant2 = pattern(7); %RH
R22(end).product1 = pattern(10); % ROOH
R22(end).product2 = pattern_non_reactive(1); % R.
R22(end).k = 0;
R22(end).type = 12;

% 22.4 Hydrogen abstraction from C_BA by OO. k_H_dc
R22(end+1).reactant1 = pattern(9); % ROO.
R22(end).reactant2 = pattern(7); %RH
R22(end).product1 = pattern(10); % ROOH
R22(end).product2 = pattern_non_reactive(2); % R.
R22(end).k = 0;
R22(end).type = 12;
 
% 22.5 Hydrogen abstraction from C_BA by O. k_H_dc
R22(end+1).reactant1 = pattern(11); % RO.
R22(end).reactant2 = pattern(7); %RH
R22(end).product1 = pattern(13); % ROH
R22(end).product2 = pattern_non_reactive(1); %R.
R22(end).k = 0;
R22(end).type = 13;


% 22.6 Hydrogen abstraction from C_BA by O. k_H_dc
R22(end+1).reactant1 = pattern(11); % RO.
R22(end).reactant2 = pattern(7); %RH
R22(end).product1 = pattern(13); % ROH
R22(end).product2 = pattern_non_reactive(2); %R.
R22(end).k = 0;
R22(end).type = 13;

% 22.7 Hydroperoxide decomposition with catalyst 
R22(end+1).reactant1 = pattern(10); % ROOH
R22(end).reactant2 = pattern(1); %Co(II)
R22(end).product1 = pattern(11); % RO.
R22(end).product2 = pattern(5); %Co(III)OH
R22(end).k = 0;
R22(end).type = 14;

% 22.8 Hydroperoxide decomposition with catalyst
R22(end+1).reactant1 = pattern(31); % R(=O)OOH
R22(end).reactant2 = pattern(1); %Co(II)
R22(end).product1 = pattern(29); % R(=O)O.
R22(end).product2 = pattern(5); %Co(III)OH
R22(end).k = 0;
R22(end).type = 14;


%__________________________________
% 22.9 Recombination: R. + R. ->  R-R (crosslink on C2)  k_c
R22(end+1).reactant1 = pattern(8); % R.
R22(end).reactant2 = pattern(8); % R.
R22(end).product1 = pattern_non_reactive(13); % R-...
R22(end).product2 = pattern_non_reactive(13); % R-...
R22(end).k = 0;
R22(end).type = 15;

% 22.10 Recombination: R. + ROO. -> R-O-O-R (crosslink on C2) k_c
R22(end+1).reactant1 = pattern(8); % R.
R22(end).reactant2 = pattern(9); % ROO.
R22(end).product1 = pattern_non_reactive(11); % R-OO-...
R22(end).product2 = pattern(14); % R-OO-...
R22(end).k = 0;
R22(end).type = 16;

% 22.11 Recombination: RO. + RO. -> R-O-O-R   k_c
R22(end+1).reactant1 = pattern(11); % RO.
R22(end).reactant2 = pattern(11); % RO.
R22(end).product1 = pattern(14); % R-OO-...
R22(end).product2 = pattern(14); % R-OO-...
R22(end).k = 0;
R22(end).type = 17;

% 22.12 Recombination: RO. + R. -> R-O-R (crosslink on C2)  k_c
R22(end+1).reactant1 = pattern(11); % RO.
R22(end).reactant2 = pattern(8); % R.
R22(end).product1 =pattern(15); % R-O-...
R22(end).product2 = pattern_non_reactive(12); % R-O-...
R22(end).k = 0;
R22(end).type = 18;



%22.13 Hydrogen abstraction from dienal by ROO.
%_______________________________________ 
R22(end+1).reactant1 = pattern(18); % R(O)H
R22(end).reactant2 = pattern(9); % ROO.
R22(end).product1 = pattern(20); % R(O).
R22(end).product2 = pattern(10); % ROOH
R22(end).k = 0;
R22(end).type = 19;

% 22.14 Hydrogen abstraction from dienal by R(O)OO.
R22(end+1).reactant1 = pattern(18); % R(O)H
R22(end).reactant2 = pattern(21); % R(O)OO.
R22(end).product1 = pattern(20); % R(O).
R22(end).product2 = pattern(31); % R(O)OOH
R22(end).k = 0;
R22(end).type = 19;

%22.15 Bis-allylic hydrogen abstraction by peroxide after beta-scission
%_____________________________________________k_H_ba
R22(end+1).reactant1 = pattern(21); % R(O)OO.
R22(end).reactant2 = pattern(7); % RH
R22(end).product1 = pattern(31); %R(O)OOH
R22(end).product2 = pattern_non_reactive(1); %R.
R22(end).k = 0;
R22(end).type = 12;

%22.16 Bis-allylic hydrogen abstraction by peroxide after beta-scission

R22(end+1).reactant1 = pattern(21); % R(O)OO.
R22(end).reactant2 = pattern(7); % RH
R22(end).product1 = pattern(31); %R(O)OOH
R22(end).product2 = pattern_non_reactive(2); %R.
R22(end).k = 0;
R22(end).type = 12;

%22.17 Bis-allylic hydrogen abstraction by alkoxy after beta-scission
R22(end+1).reactant1 = pattern(29); % R(O)O.
R22(end).reactant2 = pattern(7); % RH
R22(end).product1 = pattern(23); %R(O)OH
R22(end).product2 = pattern_non_reactive(1); %R.
R22(end).k = 0;
R22(end).type = 13;

%22.18 Bis-allylic hydrogen abstraction by alkoxy after beta-scission
R22(end+1).reactant1 = pattern(29); % R(O)O.
R22(end).reactant2 = pattern(7); % RH
R22(end).product1 = pattern(23); %R(O)OH
R22(end).product2 = pattern_non_reactive(2); %R.
R22(end).k = 0;
R22(end).type = 13;

%22.19 
R22(end+1).reactant1 = pattern(18); % R(O)H
R22(end).reactant2 = pattern(11); % RO.
R22(end).product1 = pattern(20); %R(O).
R22(end).product2 = pattern(13); %ROH
R22(end).k = 0;
R22(end).type = 20;


%22.20
R22(end+1).reactant1 = pattern(18); % R(O)H
R22(end).reactant2 = pattern(29); % R(O)O.
R22(end).product1 = pattern(20); %R(O).
R22(end).product2 = pattern(23); %R(O)OH
R22(end).k = 0;
R22(end).type = 20;



%______________________________________
%22.21 Carboxylic acid formation Bayer Villiger
R22(end+1).reactant1 = pattern(22); % R(O)OOH
R22(end).reactant2 = pattern(18); % R(O)H
R22(end).product1 = pattern(23); %R(O)OH
R22(end).product2 = pattern_non_reactive(6); %R(O)OH
R22(end).k = 0;
R22(end).type = 21;

%______________________________________________
% termination of alkyl radical after beta-scission

% Radical reactions of pattern 22

%22.22 termination of R'. with R.
R22(end+1).reactant1 = pattern(19); % R'.
R22(end).reactant2 = pattern(8); % R.
R22(end).product1 = pattern_non_reactive(8); %R'-...
R22(end).product2 = pattern_non_reactive(13); %R-...
R22(end).k = 0;
R22(end).type = 15;


%22.23 termination of R'. with ROO.
R22(end+1).reactant1 = pattern(19); % R'.
R22(end).reactant2 = pattern(9); % ROO.
R22(end).product1 = pattern_non_reactive(9); %R'-OO-...
R22(end).product2 = pattern(14); %R-OO-...
R22(end).k = 0;
R22(end).type = 16;


%22.24 termination of R'. with RO.
R22(end+1).reactant1 = pattern(19); % R'.
R22(end).reactant2 = pattern(11); % RO.
R22(end).product1 = pattern_non_reactive(14); %R'-O-...
R22(end).product2 = pattern(15); %R-O-...
R22(end).k = 0;
R22(end).type = 18;


% 22.25 termination of R'. with R'.
R22(end+1).reactant1 = pattern(19); % R'. 
R22(end).reactant2 = pattern(19); % R'.
R22(end).product1 = pattern_non_reactive(8); % 
R22(end).product2 = pattern_non_reactive(8); % 
R22(end).k = 0;
R22(end).type = 15;
 
%__________________________________new_15.05
% automatizing recombinations of beta scission products
% 22.26-28
% pattern 20 C.
recombination_C_C( 20, 26 ); % 3 reactions
% 22.29
recombination_C_O( 20, 28 ); % 2 reactions
% 22.30
recombination_C_OO(20, 27 ); % 2 reactions
% 22.31-33 pattern 21 OO.
recombination_OO_C(21, 27 ); % 3 reactions
% 22.34-36 pattern 29 O.
recombination_O_C( 29, 28 ); % 3 reactions
% 22.37-38
recombination_O_O( 29, 27 ); % 2 reactions

% 22.39 Hydroperoxides becomes cobaltperoxide
R22(end+1).reactant1 = pattern(10); % ROOH
R22(end).reactant2 = pattern(5); % Co(III)OH
R22(end).product1 = pattern(25); % ROOCo
R22(end).product2 = pattern(6); %H2O
R22(end).k = 0;
R22(end).type = 22;

% 22.40 Hydroperoxides becomes cobaltperoxide
R22(end+1).reactant1 = pattern(31); % R(=O)OOH
R22(end).reactant2 = pattern(5); %Co(III)OH
R22(end).product1 = pattern(32); % R(=O)OOCo
R22(end).product2 = pattern(6); %H2O
R22(end).k = 0;
R22(end).type = 22;



% second order rules (2 -> 3)
% __________________________________
% 23.1 Recombination: COO. + COO. -> C-O-O-C  + O2 k_c

R23(1).reactant1 = pattern(9); % ROO.
R23(1).reactant2 = pattern(9); % ROO.
R23(1).product1 = pattern(14); % R-OO-...
R23(1).product2 = pattern(14); % R-OO-...
R23(1).product3 = pattern(4); % O2
R23(1).k = 0;
R23(1).type = 26;

% 23.2 Recombination: COO. + COO. -> C-O-O-C  + O2 k_c
R23(end+1).reactant1 = pattern(21); % RC(=O)OO.
R23(end).reactant2 = pattern(9); % ROO.
R23(end).product1 = pattern(27); % RC(=O)-OO-...
R23(end).product2 = pattern(14); %R-OO-...
R23(end).product3 = pattern(4); %O2
R23(end).k = 0;
R23(end).type = 26;

% 23.3 Recombination: COO. + COO. -> C-O-O-C  + O2 k_c
R23(end+1).reactant1 = pattern(21); % RC(=O)OO.
R23(end).reactant2 = pattern(21); % RC(=O)OO.
R23(end).product1 = pattern(27); % RC(=O)-OO-...
R23(end).product2 = pattern(27); % RC(=O)-OO-...
R23(end).product3 = pattern(4); %O2
R23(end).k = 0;
R23(end).type = 26;

% 23.4 Rearangement for Hydroperoxide decomposition with catalyst
R23(end+1).reactant1 = pattern(25); % ROOCo
R23(end).reactant2 = pattern(10); % ROOH
R23(end).product1 = pattern(11); % RO.
R23(end).product2 = pattern(9); % ROO.
R23(end).product3 = pattern(5); %Co(III)OH
R23(end).k = 0;
R23(end).type = 27;

% 23.5 Rearangement for Hydroperoxide decomposition with catalyst
R23(end+1).reactant1 = pattern(25); % ROOCo
R23(end).reactant2 = pattern(31); % R(=O)OOH
R23(end).product1 = pattern(11); % RO.
R23(end).product2 = pattern(21); % R(=O)OO.
R23(end).product3 = pattern(5); %Co(III)OH
R23(end).k = 0;
R23(end).type = 27;

% 23.6 Rearangement for Hydroperoxide decomposition with catalyst
R23(end+1).reactant1 = pattern(32); % R(=O)OOCo
R23(end).reactant2 = pattern(10); % ROOH
R23(end).product1 = pattern(29); % R(=O)O.
R23(end).product2 = pattern(9); % ROO.
R23(end).product3 = pattern(5); %Co(III)OH
R23(end).k = 0;
R23(end).type = 27;


% 23.7 Rearangement for Hydroperoxide decomposition with catalyst
R23(end+1).reactant1 = pattern(32); % R(=O)OOCo
R23(end).reactant2 = pattern(31); % R(=O)OOH
R23(end).product1 = pattern(29); % R(=O)O.
R23(end).product2 = pattern(21); % R(=O)OO.
R23(end).product3 = pattern(5); %Co(III)OH
R23(end).k = 0;
R23(end).type = 27;


% % 23.8 russell termination
% 
% R23(end+1).reactant1 = pattern(30); % RHCOO.
% R23(end).reactant2 = pattern(30); % RHCOO.
% R23(end).product1 = pattern_non_reactive(21); % RC=O
% R23(end).product2 = pattern_non_reactive(22); % RHCOH
% R23(end).product3 = pattern(4); %O2
% R23(end).k = 0;
% R23(end).type = 28;

end









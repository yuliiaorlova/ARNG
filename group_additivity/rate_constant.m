function [ value ] = rate_constant( T, type, Hrxn )

R = 1.987e-3; % (cal/mole)
%T = 293;

rate_parameters = struct('A', [], 'Ea', [], 'alpha', []);

% Initiation Co(III)OOH -> Co(II) type 2
rate_parameters(1).A = 0.5;
rate_parameters(1).Ea = 0;
rate_parameters(1).alpha = 0;

% Epoxide ring opening type 3
rate_parameters(2).A = 0;
rate_parameters(2).Ea = 0;
rate_parameters(2).alpha = 0;

% Cobalt peroxide decomposition type 4
rate_parameters(3).A = 0.0095;
rate_parameters(3).Ea = 0;
rate_parameters(3).alpha = 0;

% Hydroperoxide decomposition type 5
rate_parameters(4).A = 1e15; %5e8;
rate_parameters(4).Ea = 0; %40;
rate_parameters(4).alpha = 1;

% Beta scission D1 type 6
rate_parameters(5).A = 1e14;
rate_parameters(5).Ea = 9.5;
rate_parameters(5).alpha = 0.85;

% Beta scission D2 type 7
rate_parameters(6).A = 1e14;
rate_parameters(6).Ea = 8.7;
rate_parameters(6).alpha = 0.85;

% Epoxidation from hydroperoxide type 8
rate_parameters(7).A = 1e13;
rate_parameters(7).Ea = 19.3;
rate_parameters(7).alpha = 0.8;

% Oxidation of catalyst (primary initiation) type 9
rate_parameters(8).A = 0.666; % 1e13
rate_parameters(8).Ea = 0; % heat of reaction
rate_parameters(8).alpha = 0;

% C radical oxidation type 10
rate_parameters(9).A = 1e8;
rate_parameters(9).Ea = 0;
rate_parameters(9).alpha = 0;

% Initiation with Co: CH + Co(III)OO.  (secondary initiation) 11
rate_parameters(10).A = 0.5*5e-5; % 1e15
rate_parameters(10).Ea = 0; % heat of reaction
rate_parameters(10).alpha = 0;

% Bis allylic H abstraction by OO. type 12
rate_parameters(11).A = 6.6; % 0.5*1e7; 
rate_parameters(11).Ea = 0; % 8.42
rate_parameters(11).alpha = 0;

% Bis allylic H abstraction by O. type 13
rate_parameters(12).A = 0.5*1e6;
rate_parameters(12).Ea = 11.9;
rate_parameters(12).alpha = 0.91;

% Hydroperoxide decomposition with catalyst 14
rate_parameters(13).A = 2.7e8; % for acac, 1.7e14 for oct
rate_parameters(13).Ea = 13; % for acac, 19.7 for oct
rate_parameters(13).alpha = 0;

% Recombination C. + C. type 15
rate_parameters(14).A = 1e8;
rate_parameters(14).Ea = R*T;
rate_parameters(14).alpha = 0;

% Recombination C. + COO. 16
rate_parameters(15).A = 1e8;
rate_parameters(15).Ea = R*T;
rate_parameters(15).alpha = 0;

% Recombination CO. + CO. type 17
rate_parameters(16).A = 1e8;
rate_parameters(16).Ea = R*T;
rate_parameters(16).alpha = 0;

% Recombination CO. + C. type 18
rate_parameters(17).A = 1e8;
rate_parameters(17).Ea = R*T;
rate_parameters(17).alpha = 0;

% Hydrogen abstraction from dienal R(O)H by ROO. type 19
rate_parameters(18).A = 1e7;
rate_parameters(18).Ea = 6.2*2;
rate_parameters(18).alpha = 1;

% Hydrogen abstraction from dienal R(O)H by RO. 20
rate_parameters(19).A = 1e7;
rate_parameters(19).Ea = 11.9;
rate_parameters(19).alpha = 0.91;

% Bayer Villiger type 21
rate_parameters(20).A = 1.2e4;
rate_parameters(20).Ea = 8.5;
rate_parameters(20).alpha = 0;

% Hydroperoxides becomes cobaltperoxide   type 22
rate_parameters(21).A = 0.0635;
rate_parameters(21).Ea = 0;
rate_parameters(21).alpha = 0;

% Peroxy crosslink decomposition type 23
rate_parameters(22).A = 1e15; %5e8;
rate_parameters(22).Ea = 0; %40;
rate_parameters(22).alpha = 1;

% Addition type 24
rate_parameters(23).A = 1e8; % 0.5702; 
rate_parameters(23).Ea = 11.24;
rate_parameters(23).alpha = 0.24;

% Epoxidation from R-OO-R. type 25
rate_parameters(24).A = 1e11;
rate_parameters(24).Ea = 17.8;
rate_parameters(24).alpha = 0.8;

% Recombination ROO. + ROO. type 26
rate_parameters(25).A = 1e6;
rate_parameters(25).Ea = R*T;
rate_parameters(25).alpha = 0;

% Rearrangement for Hydroperoxide decomposition with catalyst type 27
rate_parameters(26).A = 1.7e14*1.7e14; % (1st step * 2nd step)
rate_parameters(26).Ea = 22.22+5.9; % (1st step + 2nd step)
rate_parameters(26).alpha = 0;

% Russell termination type 28
rate_parameters(27).A = 1e5;
rate_parameters(27).Ea = 0;
rate_parameters(27).alpha = 0;

% Hydrogen transfer type 29
rate_parameters(28).A = 1e7;
rate_parameters(28).Ea = 9.1;
rate_parameters(28).alpha = 0.7;

if rate_parameters(type-1).Ea == 0 && rate_parameters(type-1).alpha == 0
    value = rate_parameters(type-1).A;
else
    
%    if Hrxn<0
        if (rate_parameters(type-1).Ea+rate_parameters(type-1).alpha*Hrxn) <0
            Hrxn = 0;
            value = rate_parameters(type-1).A*exp((-(rate_parameters(type-1).Ea+rate_parameters(type-1).alpha*Hrxn))/(R*T));
        else
            value = rate_parameters(type-1).A*exp((-(rate_parameters(type-1).Ea+rate_parameters(type-1).alpha*Hrxn))/(R*T));
        end
%      else
%         if (rate_parameters(type-1).Ea+(1-rate_parameters(type-1).alpha)*Hrxn) < 0
%             Hrxn = 0;
%             value = rate_parameters(type-1).A*exp((-(rate_parameters(type-1).Ea+(1-rate_parameters(type-1).alpha)*Hrxn))/(R*T));
%         else
%             value = rate_parameters(type-1).A*exp((-(rate_parameters(type-1).Ea+(1-rate_parameters(type-1).alpha)*Hrxn))/(R*T));
%         end
%     end
    
    if value<0
        value = 0;
    end
end
end


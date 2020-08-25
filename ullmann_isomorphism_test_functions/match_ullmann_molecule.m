function [ idx ] = match_ullmann_molecule( A, labels, stuff )
%MATCH_ULLMANN_MOLECULE check if molecule is already in the list, if so -
%returns its index in the list
global Species;
idx = [];
indicator = 0;
if length(A)<3
    for i = 1 : length(Species)
        if length(Species(i).adj)<3
            if isequal(Species(i).lbl, labels)
                idx = i;
                break
            end
        end
    end
else
    %     % some basic statistics on A in order to cut off some cases for comparison
    %     n1 = length(A);
    %     O1 = sum(ismember(labels, 'O'));
    %     H1 = sum(ismember(labels, 'H'));
    %     C1 = sum(ismember(labels, 'C'));
    stuff = get_stuff(A, labels);
    
    
    for i = 1 : length(Species)
        
        if length(Species(i).adj)>=3
            %             n2 = length(Species(i).adj);
            %             O2 = sum(ismember(Species(i).lbl, 'O'));
            %             H2 = sum(ismember(Species(i).lbl, 'H'));
            %             C2 = sum(ismember(Species(i).lbl, 'C'));
            
            if isequal(stuff, Species(i).stuff) %n1==n2 && O1==O2 && H1==H2 && C1==C2
                [ indicator ] = ullmann_molecule( A, labels, Species(i).adj, Species(i).lbl );
            end
        end
        if indicator
            idx = i;
            break
        end
        
    end
end
end


function   graph2mol( A, labels, id)
% graph2mol converts molecular graphs into mol files

fileID = fopen(sprintf('species_%d.mol', id),'w', 'native', 'UTF-8');

% first 3 lines: 1) name of the molecule, 2) software it was created with,
% 3) comment line
fprintf(fileID,'\n');
fprintf(fileID,'\n');
fprintf(fileID,'\n');
% 4) number of atoms, number of bonds
n = length(A);
m = nnz(tril(A))-nnz(diag(A));

fprintf(fileID, '%3d%3d  0  0  0  0  0  0  0  0999 V2000\n', n, m);
for i = 1 : n
fprintf(fileID, '    0.0000    0.0000    0.0000 %s   0  0  0  0  0  0  0  0  0  0  0  0\n', labels{i});
end

S = tril(A);
S(logical(eye(size(S)))) = 0;
[row, col] = find(S>0);
for i =1  : length(row)
fprintf( fileID, '%3d%3d%3d  0  0  0\n', col(i), row(i), floor(S(row(i), col(i)))) ;
end
radicals = find(diag(A)==1);
if ~isempty(radicals)
radicals_string = sprintf('M  RAD%3d', length(radicals));
for i = 1 : length(radicals)
    
    radicals_string = strcat(radicals_string, sprintf(' %3d %3d', radicals(i), 2));

end
radicals_string = strcat(radicals_string, '\n');
fprintf(fileID,radicals_string);
end
fprintf(fileID,'M  END\n');
fclose(fileID);


end


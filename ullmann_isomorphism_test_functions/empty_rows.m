function [ indicator ] = empty_rows( M )
indicator = 0;
s = sum(M, 2);
if ~isempty(find(s==0))
    indicator = 1;
end

end


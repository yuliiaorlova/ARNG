function stuff = get_stuff(A, labels)
    stuff = [];
    if length(A)>1
    spec=round(eigs(balance(A),length(A)) * 1e5 );
    spec = sort(spec);
    labels =sort(double(cell2mat(labels')));
    stuff = [spec' labels ];
    end
   
end

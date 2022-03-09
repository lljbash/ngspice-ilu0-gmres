mkdir('ref_sol');
for i = 0 : 99
    filename = sprintf('MM1/%02d', i);
    tmp = load(filename);
    n = tmp(1, 1);
    tmp = [tmp(2:end, :); n n 0];
    A = spconvert(tmp);
    filename = sprintf('rhs/%02d', i);
    rhs = load(filename);
    x = A \ rhs;
    filename = sprintf('ref_sol/%02d', i);
    dlmwrite(filename, x);
end

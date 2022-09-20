function hxMat = generateHistoryMatrix(input, nBack)

% transpose to be a column vector
if size(input, 1) == 1
    input = input';
end

hxMat = NaN(length(input), nBack);
for t = 1:nBack
    hxMat(:, t) = [NaN(t, 1); input(1:end - t)];
end
function y = cell_times(t,lda)
% CELL_TIMES(t,lda) scalar t times matrix array lda
y = cell(size(lda));
for i = 1:size(lda,1)
    for j = 1:size(lda,2)
        y{i,j} = lda{i,j} * t;
    end
end
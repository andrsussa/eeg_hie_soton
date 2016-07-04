function ploteegrel( M )

    n = size(M,1);
    L = num2cell(1:n);
    imagesc(M(:,:,1));
    set(gca, 'XTick', 1:n);
    set(gca, 'YTick', 1:n);
    set(gca, 'XTickLabel', L);
    set(gca, 'YTickLabel', L);
    title('COH', 'FontSize', 14);
    colormap('jet');
    colorbar;

end
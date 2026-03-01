function showedges(image, xEdges, yEdges)

hold off;
imagesc(image)
hold on;
% regs (c+0.5,r-0.5) : (c+0.5,r+0.5)
% ondr (c-0.5,r+0.5) : (c+0.5,r+0.5)
xEdges = xEdges/max(max(xEdges));
yEdges = yEdges/max(max(yEdges));
xEdges(find(xEdges < 0.5) ) = 0.0;
yEdges(find(yEdges < 0.5) ) = 0.0;

% regs
[nr,nc] = size(xEdges);
for i=1:nr
  r = i+2;
  for j = 1:nc
    c = j+2;
    if (xEdges(i,j) > 0)
      plot([c+0.4, c+0.4],[r-0.5, r+0.5], 'w');
      plot([c+0.5, c+0.5],[r-0.5, r+0.5], 'k');
      plot([c+0.6, c+0.6],[r-0.5, r+0.5], 'w');
    endif
  endfor
endfor

% onder
[nr,nc] = size(yEdges);
for i=1:nr
  r = i+2;
  for j = 1:nc
    c = j+2;
    if (yEdges(i,j) > 0)
      plot([c-0.5, c+0.5],[r+0.4, r+0.4], 'w');
      plot([c-0.5, c+0.5],[r+0.5, r+0.5], 'k');
      plot([c-0.5, c+0.5],[r+0.6, r+0.6], 'w');
    endif
  endfor
endfor

hold off;

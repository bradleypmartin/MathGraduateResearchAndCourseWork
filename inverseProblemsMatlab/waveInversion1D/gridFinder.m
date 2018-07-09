intFinder = 0;
queryPoint = 0.545001;

for m = 2:N
   if xgrid(m-1,1) < queryPoint && xgrid(m,1) > queryPoint
      intFinder = m; 
   end
end

intFinder
compare = [];
for i = 1:length(vectors)
  comparison =  [radtodeg(acos(dot(vectors(i,:),normalsJan(i,:))/norm(vectors(i,:))*norm(normalsJan(i,:))))] 
  compare = [compare; comparison]
end

tbl = array2table(compare)
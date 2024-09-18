factors=3.*2.^(0:3);
factors = [2,3,6,7,8,11,12]
nodeNumbers=zeros(size(factors));
numFaces=zeros(size(factors));
numDOF=zeros(size(factors));
for i=1:size(factors,2)
  factor=factors(i);
 [ node_num, edge_num, triangle_num ] = sphere_grid_icos_size ( factor );
 nodeNumbers(i)=node_num;
 numFaces(i)=triangle_num; 
 numDOF(i)=node_num*3;
end
% pt1=[0,1,1]'; 
% pt2=[1,1/2,0]'; 
% pt3=[1/5,1,1]';

pt1=[-1/2,0,0]';
pt2=[1/2,0,0]'; 
pt3=[0,sqrt(3)/2,0]';

xcoords=[pt1(1),pt2(1),pt3(1),pt1(1)]'; 
ycoords=[pt1(2),pt2(2),pt3(2),pt1(2)]'; 
zcoords=[pt1(3),pt2(3),pt3(3),pt1(3)]'; 

plot3(xcoords,ycoords,zcoords)

vertices=[pt1,pt2,pt3]
directions=[pt1-pt2,pt2-pt3,pt3-pt1];
lengths=vecnorm(directions); 
directions=directions./lengths;
normaltoplane=cross(directions(:,1),directions(:,2)); 
normaltoplane=normaltoplane./vecnorm(normaltoplane);
normalstosides=[cross(normaltoplane,directions(:,1)), ...
    cross(normaltoplane,directions(:,2)),...
    cross(normaltoplane,directions(:,3))]
ht1=lengths(2).*dot(directions(:,2),normalstosides(:,1));
bh=lengths(1).*ht1; 
heights=[ht1,bh./lengths(2),bh./lengths(3)];

Triangle=struct('vertices',vertices,'directions',directions, ...
    'normaltoplane',normaltoplane,'normalstosides',normalstosides,...
    'lengths',lengths,'heights',heights,'bh',bh);

disp('check areas')
disp(['base times height = ' num2str(Triangle.bh)]) 
disp(['cross prod1 = ' num2str(Triangle.lengths(1).*Triangle.lengths(2).*vecnorm(cross( ...
    Triangle.directions(:,1),Triangle.directions(:,2))))])
disp(['cross prod2 = ' num2str(Triangle.lengths(2).*Triangle.lengths(3).*vecnorm(cross( ...
    Triangle.directions(:,2),Triangle.directions(:,3))))])
disp(['cross prod3 = ' num2str(Triangle.lengths(3).*Triangle.lengths(1).*vecnorm(cross( ...
    Triangle.directions(:,3),Triangle.directions(:,1))))])

% plot3([pt1,pt2,pt3,pt1],'o-')
gd = [3, 4,x(:)',y(:)']'; % rectangle
ns = uint16('R1')';
sf = 'R1';

dl = decsg(gd,sf,ns);

model = createpde();
geometryFromEdges(model,dl);
msh = generateMesh(model,'Hmax',h,'GeometricOrder','linear');

points = msh.Nodes;
elements = msh.Elements'-1;

tr = triangulation(elements+1,points');
edges = tr.edges-1;
function intersect = TriangleSegmentIntercetion(p1,p2,v0,v1,v2)
% Check if triangle and segment interact
% Based on Triangle/Ray Intersection version 1.7 (309 KB) by Jaroslaw Tuszynski

zero = -eps;
dir = p2-p1;
intersect = 0;
edge1 = v1-v0;          % find vectors for two edges sharing vert0
edge2 = v2-v0;
tvec  = p1 - v0;          % vector from vert0 to ray origin
pvec  = cross(dir,edge2);  % begin calculating determinant - also used to calculate U parameter
det   = sum(edge1.*pvec);   % determinant of the matrix M = dot(edge1,pvec)

angleOK = (abs(det)>eps);

if all(~angleOK)
    return;
end % if all parallel than no intersections

det(~angleOK) = nan;              % change to avoid division by zero
u    = sum(tvec.*pvec)./det;    % 1st barycentric coordinate

  v = nan+zeros(size(u));
  t=v;
  ok = (angleOK & u>=-zero & u<=1.0+zero); % mask
  % if all line/plane intersections are outside the triangle than no intersections
  if ~any(ok), intersect = ok; return; end
  qvec = cross(tvec, edge1,2); % prepare to test V parameter
  v(ok,:) = sum(dir.*qvec,2) ./ det; % 2nd barycentric coordinate
  t(ok,:) = sum(edge2.*qvec,2)./det;
  % test if line/plane intersection is within the triangle
  ok = (ok & v>=-zero & u+v<=1.0+zero);

  intersect = (ok & t>=-zero & t<=1.0+zero);

end
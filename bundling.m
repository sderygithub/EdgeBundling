function comp = computeCompatiblity(edges,v)
   % Indices
   first = 1:3;
   last = 4:6;
   
   % Self implemented is faster
   %  pair = choosenk(length(ep),2);
   % Matlab version is too slow
   %  pair = nchoosek(1:length(ep),2);
   
   mid = v(edges(:,1),:) + 0.5 * normalizeVector(v(edges(:,2),:) - v(edges(:,1),:));

   comp = zeros(length(edges),length(edges));
   for iEdge = 1:length(edges)
      % 
      cEdge = ones(length(edges),1);
      cEdge(iEdge) = 0;
      if all(cEdge == 0)
         continue;
      end
      
      ei = [v(edges(iEdge,1),:), v(edges(iEdge,2),:)];
      ej = [v(edges(cEdge == 1,1),:), v(edges(cEdge == 1,2),:)];
      
      % Length
      ei_length = sqrt(sum((ei(:,first) - ei(:,last)).^2,2));
      ej_length = sqrt(sum((ej(:,first) - ej(:,last)).^2,2));

      % Midlines
      ei_mid = mid(iEdge,:);
      ej_mid = mid(cEdge == 1,:);

      % Angle compatiblity
      angle = cross(repmat((ei(:,first)-ei(:,last)),[],size(ej,1)),(ej(:,last)-ej(:,first)),2);
      angle(:,1) = angle(:,1) ./ (ei_length * ej_length);
      angle(:,2) = angle(:,2) ./ (ei_length * ej_length);
      angle(:,3) = angle(:,3) ./ (ei_length * ej_length);
      ca = 1 - sqrt(sum(angle.^2,2));
      % ca = exp(sqrt(sum(angle.^2,2))).^2;
      % ca = 1 - (ca - min(ca)) / (max(ca) - min(ca));
      
      % Scale compatibility
      lavg = (ei_length + ej_length) / 2;
      cs = 2 ./ (lavg ./ min(ei_length,ej_length) + max(ei_length,ej_length) ./ lavg);
      
      % Position compatibility
      cp = lavg ./ (lavg + sqrt(sum((ej_mid - repmat(ei_mid,[],size(ej,1))).^2,2)));
      
      % Edge compatibility
      if any(cs .* cp > 0.9)
         % Add a visibility criteria
         comp(iEdge,cEdge == 1) = ca .* cs .* cp;
      else
         comp(iEdge,cEdge == 1) = ca .* cs .* cp;
      end
   end
end

% 0-1-2-3-4-5-6-7
% 2-3-4-5-6-7-8-9
function points = bundleTracks(edges,v,points)
   % 
   nCycle = 6;
   nSmoothing = 1;
   % 
   comp = computeCompatiblity(edges,v);
   % 
   for iCycle=1:nCycle
      fprintf(sprintf('BST> Cycle %d\n',iCycle));
      np = size(points,2) / 3;
      %points = smooth(points, np - 2 + 1);
      points = subdivide(points);
      % Inverse relation to reduce attraction at later stage
      nIter = nCycle - iCycle + 1;
      for iter=1:nIter
         fprintf(sprintf('BST>   Iteration %d\n',iter));
         points = attract(points, comp, nIter);
      end
   end
   % Smoothing
   if nSmoothing > 0
      fprintf(sprintf('BST>   Smoothing operation\n'));
      for iCycle=1:nSmoothing
         points = subdivide(points);
      end
      points = attract(points, comp, nIter);
   end
end

function points = subdivide(p)
   np = size(p,2) / 3;
   points = zeros(size(p,1),(np + 1) * 3);
   points(:,1:3) = p(:,1:3);
   for i=1:np-1
      p1 = p(:,((i-1)*3+1):((i-1)*3)+3);
      p2 = p(:,((i)*3+1):((i)*3)+3);
      d = p2 - p1;
      bv = normalizeVector(d);
      bl = sqrt(sum(d.^2,2));
      points(:,(i*3+1):(i*3+3)) = p1 + bsxfun(@times, bv, 0.5 * bl);
   end
   points(:,end-2:end) = p(:,end-2:end);
end

function points = smooth(p,newp)
   np = size(p,2) / 3;
   points = zeros(size(p,1),(newp+2)*3);
   % Compute current total length
   polyLength = 0;
   for i=1:np-1
      p1 = p(:,((i-1)*3+1):((i-1)*3)+3);
      p2 = p(:,((i)*3+1):((i)*3)+3);
      polyLength = polyLength + sqrt(sum((p2 - p1).^2,2));
   end
   % New segment length
   newl = polyLength / (newp + 1);
   lengthSoFar = zeros(size(polyLength));
   % Start/End points
   points(:,1:3) = p(:,1:3);
   points(:,end-2:end) = p(:,end-2:end);
   % Compute new points
   p1 = p(:,((1-1)*3+1):((1-1)*3)+3);
   p2 = p(:,((1)*3+1):((1)*3)+3);
   seg = 1;
   for i=1:newp
      % Find before and after points
      while (segLength(p,seg) < newl*i-lengthSoFar)
         lengthSoFar = lengthSoFar + segLength(p,seg);
         seg = seg + 1;
         p1 = p(:,((seg-1)*3+1):((seg-1)*3)+3);
         p2 = p(:,((seg)*3+1):((seg)*3)+3);
      end
      ib = newl * i - lengthSoFar;
      d = p2 - p1;
      bv = normalizeVector(d);
      points(:,(i*3+1):(i*3+3)) = p1 + bsxfun(@times, bv, ib);
   end
end

function l = segLength(p, n)
   p1 = p(:,((n-1)*3+1):((n-1)*3)+3);
   p2 = p(:,((n)*3+1):((n)*3)+3);
   l = sqrt(sum((p2 - p1).^2,2));
end

function v = normalizeVector(p)
   pl = sqrt(sum(p.^2,2));
   valid = ~(isnan(pl) | pl == 0);
   v = zeros(size(p));
   v(valid == 1,1) = bsxfun(@rdivide, p(valid == 1,1), pl(valid == 1));
   v(valid == 1,2) = bsxfun(@rdivide, p(valid == 1,2), pl(valid == 1));
   v(valid == 1,3) = bsxfun(@rdivide, p(valid == 1,3), pl(valid == 1));
end

function p = attract(p,comp,w)
   if nargin < 3
      w = 2;
   end
   % Constant
   c_thr = 0.5;
   % Compute forces
   np = size(p,2) / 3;
   force = zeros(size(p));
   % For all edges
   for ie=1:size(p,1)
      % For every point (except start/end)
      for ip=2:np-1
         % Take first vector
         pi = p(ie,((ip-1)*3+1):(ip*3));
         % Retrieve attractive vectors
         pe = p((comp(ie,:) > c_thr),((ip-1)*3+1):(ip*3));
         % Duplicate pi for convenience
         pempi = pe - repmat(pi,size(pe,1),1);
         % Distance
         de = sqrt(sum(pempi.^2,2));
         % Weight definition
         weight = normpdf(de, 4, 2) * w;
         % Prevent "/ by 0"
         if (size(weight,1) > 0)
            % Mean of the attractive force
            force(ie,((ip-1)*3+1):(ip*3)) = sum(bsxfun(@times,normalizeVector(pempi),weight),1) / size(weight,1);
         end
      end
   end
   % Apply force
   p = p + force;
end
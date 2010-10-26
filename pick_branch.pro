function pick_branch, x, y, xloc, height, clusters
  nst = n_elements(height)
  nleaf = n_elements(clusters[0,*]) + 1
  hi = height
  lo = hi
  lo[clusters]=rebin(1#height[nleaf:*], 2, nleaf-1)
  
  candidate = (y gt lo) and (y lt hi)

  dist = abs(x - xloc)
  dist += max(dist)*(~candidate)

  best = min(dist, result)
  return, result
end

function dendrocpp2idl, file, levels = levels
  if n_params() ne 1 then begin
     print, 'calling sequence:'
     print, 'result = dendrocpp2idl(file)'
     return, !vaues.f_nan
  endif

  if size(file, /type) ne 7 then $
     message, 'file must be a string'

  if ~file_test(file) then $
     message, 'file not found: '+file

  im =mrdfits(file, 0, h,/silent)
  id =mrdfits(file, 1, h,/silent)
  clusters = mrdfits(file, 2, h,/silent)
  seeds = reform(mrdfits(file, 3, h, /silent))
  assert, array_equal(seeds[sort(seeds)], seeds), 'Seeds are not in sorted order!'

  cv = dendrocpp2cloudviz(file)
  nseed = n_elements(seeds)

  ci = lindgen(n_elements((*cv).value))
  sz = size(im)
  ind = lindgen(n_elements(im))
  xx = ind mod sz[1]
  yy = (ind / sz[1]) mod sz[2]
  zz = (sz[0] eq 3) ? (ind / (sz[1] * sz[2])) : ind*0
  order = sort((*cv).xlocation[0:nseed-1])
  if n_elements(levels) eq 0 then begin
     levels = (*cv).height
     levels = levels[uniq(levels, sort(levels))]
     doflat = 0
  endif else doflat = 1

  ;- make the merger matrix
  merger = fltarr(nseed, nseed) + min((*cv).height)
  height = (*cv).height
  clusters = (*cv).clusters
  ind = indgen(nseed)
  merger[ind, ind] = height[0:nseed-1]
  for i = 0, n_elements(clusters[0,*]) - 1, 1 do begin
     left = clusters[0, i]
     right = clusters[1,i]
     l1 = leafward_mergers(left, clusters)
     l2 = leafward_mergers(right, clusters)
     l1 = l1[where(l1 lt nseed)]
     l2 = l2[where(l2 lt nseed)]
     for j = 0, n_elements(l1) - 1, 1 do begin
        merger[l1[j], l2] >= height[i + nseed]
        merger[l2, l1[j]] >= height[i + nseed]
     endfor
  endfor

  ;- newmerger's rows and columns are
  ;- transposed so that kernel[order[i]] maps to row/column i
  newmerger = merger
  indices, newmerger, x, y
  newmerger = newmerger[rebin(order, nseed, nseed), y]
  newmerger = newmerger[x, rebin(1#order, nseed, nseed)]

  v = (*cv).value
  l = (*cv).cluster_label
  if doflat then begin
     v = reform(v, n_elements(v), /over)
     l = reform(l, n_elements(l), /over)
  endif

  st = {$
       t: v, $
       clusters: (*cv).clusters, $
       cluster_label: l, $
       cluster_label_h: (*cv).cluster_label_h, $
       cluster_label_ri: (*cv).cluster_label_ri, $
       xlocation: (*cv).xlocation, $
       height: (*cv).height, $
       cubeindex: ci, $
       x:xx, y:yy, v:zz, $
       szdata: sz, $
       sz: sz, $
       seeds: seeds, $
       order: order, $
       err: v * 0, $
       levels: levels, $
       kernels: seeds, $
       merger: merger, $
       newmerger: newmerger $
       }

  ptr_free, cv
  return, ptr_new(st, /no_copy)
end

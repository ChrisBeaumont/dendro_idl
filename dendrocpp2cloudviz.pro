function dendrocpp2cloudviz, file

  if n_params() ne 1 then begin
     print, 'calling sequence:'
     print, ' result = dendrocpp2cloudviz(file)'
     return, !values.f_nan
  endif

  if ~file_test(file) then $
     message, 'File not found: ' + file

  im =mrdfits(file, 0, h,/silent)
  id =mrdfits(file, 1, h,/silent)
  clusters = mrdfits(file, 2, h,/silent)
  seeds = reform(mrdfits(file, 3, h,/silent))
  assert, array_equal(seeds[sort(seeds)], seeds), 'Seeds are not in sorted order!'

  sz = size(clusters)
  start = (sz[2]+1)/2
  
  assert, max(clusters[*, start-1]) eq -1 && min(clusters[*,start]) ge 0
  clusters = clusters[*, (sz[2]+1) / 2 : *]
  nleaf = (sz[2]+1)/2
  heights = replicate(-1 * !values.f_infinity, sz[2])
  xlocation = heights * 0

  h = histogram(id, min = 0, rev = ri, max = 2 * nleaf - 1)

  st = {$
       value: im, $
       clusters: clusters, $
       cluster_label:id, $
       cluster_label_h:h, $
       cluster_label_ri:ri, $
       xlocation:xlocation, $
       height:heights, $
       seeds: seeds}
  
  ptr = ptr_new(st, /no_copy)

  index = obj_new('dendro_index', ptr)

  ;- compute brightest/faintest pixels for each histogram bin
  bright = replicate(!values.f_nan, n_elements(h))
  faint = replicate(!values.f_nan, n_elements(h))
  for i = 0, n_elements(h) - 1, 1 do begin
     if h[i] eq 0 then continue
     ind = ri[ri[i] : ri[i+1] - 1]
     bright[i] = max(im[ind], /nan, min = lo)
     faint[i] = lo
  endfor

  for i = 0, sz[2] - 1, 1 do begin
     isLeaf = i lt nleaf

     ;- height of leaf is brightest pixel in that leaf
     if isLeaf then begin
        assert, h[i] ne 0
        assert, finite(bright[i])
        (*ptr).height[i] = bright[i]
     ;- height of branch is dimmest pixel in either sub-branch
     endif else begin
        left = index->leafward_mergers(clusters[0, i-nleaf])
        right = index->leafward_mergers(clusters[1, i-nleaf])

        (*ptr).height[i] = min([faint[left], faint[right]], /nan)
        assert, finite( (*ptr).height[i])
        assert, (*ptr).height[i] le min( (*ptr).height[[left,right]]), $
                'Bad tree -- merger height above children'
     endelse
  endfor

  linkDistance = max((*ptr).height) - (*ptr).height
  
  dendrogram_mod, clusters, linkDistance, ov, oc, xlocation = xlocation
  (*ptr).xlocation = xlocation

  obj_destroy, index
  return, ptr
end
       
pro test
  ptr = dendrocpp2cloudviz('~/Documents/workspace/dendro/DEBUG/ppv_big_dendro.fits')
  (*ptr).height = alog10((*ptr).height)
  dendroviz, ptr
end

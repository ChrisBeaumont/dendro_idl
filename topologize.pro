pro topologize, data, mask, $
                friends = friends, $
                specfriends = specfriends, $
                all_neighbors = all_neighbors, $
                kernels = kernels, $
                delta = delta, $
                minpix = minpix, $
                minpeak = minpeak, $
                fast = fast, $
                contour_res = contour_res, $
                levels = levels, $
                nlevels = nlevels, $
                pointer = pointer, $
                structure = structure
;+
; NAME:
;    TOPOLOGIZE
; PURPOSE:
;    Topologize evaulates the evolving connectivity of a 2 or 3 dimensional
;    data set as it is contoured on a range of levels.  It is used to
;    evaluate the heirarchical structure in an image / data cube.
;
; CALLING SEQUENCE:
;    TOPOLOGIZE, data, mask
; INPUTS:
;    DATA -- a 2 or 3 dimensional image.
;    MASK -- A binary mask array of the same dimensions as DATA that
;            indicates where the data should be studied (1=include,
;            0=exclude). 
;
; KEYWORD PARAMETERS:
;    FRIENDS, SPECFRIENDS -- Used to establish the local maxima in the
;                            data cube if not specified.  See
;                            ALLLOCMAX.pro for details. Defaults to 3.
;    ALL_NEIGHBORS -- Whether to consider corner pixels as connected to a
;                     given pixel. Defaults to zero.
;    LEVELS -- Input/output keyword.  Sepecified the contouring levels
;              at which to evaluate the data. If /FAST is set and
;              Levels is not provided, then the variable will hold the
;              ordered list of merger contours. If /FAST is set and
;              Levels are provided, then the dendrogram will be the
;              same, but subsequent routines (e.g. levelprops)
;              dependent on the value of levels will run differently.
;    NLEVELS -- The number of contour levels to use. Only used if /FAST is not
;               set, and LEVELS is not provided.
;    MINPIX -- Used in culling local maxima (see DECIMATE_KERNELS.pro
;              for details). Defaults to 4. Currently is not implemented when
;              /FAST is set.
;    DETLA --  Used in culling local maxima (see DECIMATE_KERNELS.pro
;              for details). Must be provided.
;    FAST -- If set, use some different procedures to speed up
;            execution time.
;    CONTOUR_RES -- The (absolute) precision to which determine the contour
;                   levels at which substructures split/merge. Only used if
;                   /FAST is set.
;
; OUTPUTS: Keyword controlled.  IF YOU WANT ANY OUTPUT, you must set a
; keyword.  You may set either:
;    STRUCTURE -- A structure containing information about the
;                 topology of the data.  Used by all other routines. See
;                 README file for description of tags
;    POINTER -- Named pointer to a heap variable containing all the
;               information that is in structure.  This is the
;               preferred method of dealing with data since it
;               accommodates multiple clouds in the data.
;
; SIDE EFFECTS:  HEAP VARIABLE LEAKAGE IF USED IMPROPERLY!
;
; MODIFICATION HISTORY:
;
;	Fri Oct 20 18:52:58 2006, Erik
;       Feb 2010: Added /FAST keyword. Chris Beaumont 
;       June 2010: Added histogram tags to the output structure. cnb
;       June 24 2010: Removed resolvetop, error, ecube keywords. added
;                     contour_res keyword. cnb.
;       Sep 9 2010: Default levels variable now holds the merger contour
;                   levels when called with /FAST. Gets along better
;                   with levelprops. cnb.
;       October 2010: bugfix. cluster_label_h now has a bin for every
;                     dendrogram id. cnb.
;-

  compile_opt idl2
  szdata = size(data)

  ;- check for proper input
  if n_elements(friends) eq 0 then friends = 3
  if n_elements(specfriends) eq 0 then specfriends = 3
  if szdata[0] eq 2 then specfriends = 0
  if not(keyword_set(all_neighbors)) then  all_neighbors = 0b
  if n_elements(minpix) eq 0 then minpix = 4
  if n_elements(delta) eq 0 then message, 'Please provide a value for delta'
  if keyword_set(fast) && n_elements(contour_res) eq 0 then message, 'Please provide a value for contour_res'
  if keyword_set(fast) && contour_res gt delta then $
     message, /continue, 'WARNING: contour resolution is greater than delta. This could lead to weird behavior...'

  if ~arg_present(structure) && ~arg_present(pointer) then $
     message, 'No output variables are provided!'
  
  ; Begin with data file and mask and reduce to minimum working size.
  ;- note, masked data should contain no nans
  if n_elements(mask) eq 0 then begin
     if total(~finite(data)) ne 0 then message, 'All nans in input must be masked'
  endif else begin
     if total(~finite(data[where(mask)])) ne 0 then message, 'All nans in input must be masked'
  endelse

  vectorify, data, mask = mask, x = x, y = y, v = v, t = t, $
             ind = cubeindex

  ; Trim to minimum size  
  cubify, x, y, v, t, cube = minicube, pad = (friends > specfriends), $
          twod = (szdata[0] eq 2), indvec = cubeindex, indcube = indcube
  if n_elements(kernels) gt 0 then begin
     newkern = kernels*0
     for i = 0, n_elements(kernels)-1 do newkern[i] = where(indcube eq kernels[i])
  endif 

  ; Establish Contouring levels if unset
  explicitLevels = n_elements(levels) ne 0
  if n_elements(levels) eq 0 then begin
     if n_elements(nlevels) eq 0 then  nlevels = (n_elements(x)/50 > 250) < 500
     levels = (max(t)-min(t))/(nlevels)*(findgen(nlevels))+min(t)
     levels = [0.0, levels]
  endif

  ; Establish local maxima if unset
  if n_elements(newkern) eq 0 then begin 
     if keyword_set(fast) then begin
        lmax = cnb_alllocmax(minicube, friends = friends, $
                             specfriends = specfriends)
        if n_elements(minpeak) ne 0 then begin
           good = where(minicube[lmax] gt minpeak, ct)
           if ct eq 0 then $
              message, 'No local maxima satisfied criteria'
           lmax = lmax[good]
        endif
        kernels = cnb_decimate_kernels(lmax, minicube, $
                                       all_neighbors = all_neighbors $
                                       , delta = delta, sigma = 1.0)
     endif else begin
        lmax = alllocmax(minicube, friends = friends, $
                         specfriends = specfriends)
        help, lmax
        if n_elements(minpeak) ne 0 then begin
           good = where(minicube[lmax] gt minpeak, ct)
           if ct eq 0 then $
              message, 'No local maxima satisfied criteria'
           lmax = lmax[good]
        endif
        
        kernels = decimate_kernels(lmax, minicube, $
                                   all_neighbors = all_neighbors $
                                   , delta = delta, sigma = 1.0 $
                                   , minpix = minpix)
     endelse 
  endif else kernels = newkern
  message, 'Kernels used in decomposition: '+$
           string(n_elements(kernels)), /con

  ; Calculate toplogy of the data cube
  if keyword_set(fast) then begin
     merger = cnb_mergefind(minicube, kernels, $
                            all_neighbors = all_neighbors, $
                            contour_res = contour_res, npix = npix)
  endif else begin
     merger = mergefind(minicube, kernels, levels = levels, $
                        all_neighbors = all_neighbors)
  endelse 
  nk = n_elements(kernels)
  assert, max(abs(minicube[kernels] - merger[indgen(nk), indgen(nk)])) lt 1e-4

  disconnected = where(merger ne merger, discct)
  if discct gt 0 then merger[disconnected] = 0.0

  ;- eliminate insignificant kernels
  if keyword_set(fast) && n_elements(newkern) eq 0 then begin
     decimate_merger, merger, kernels, minicube, npix, $
                      delta = delta, sigma = 1.0, $
                      minpix = minpix
  endif
  nk = n_elements(kernels)
  assert, max(abs(minicube[kernels] - merger[indgen(nk), indgen(nk)])) lt 1e-4

  ;- update levels variable if /FAST is set
  if keyword_set(fast) && ~explicitLevels then $
     levels = merger[uniq(merger, sort(merger))]

  ; Turn into sparse values again.
  ;- note that indcube[newx, newy, newv] still maps onto cubeindex, 
  ;- which itself maps into the original data. The following are true, 
  ;- even after the next line of code is executed:
  ;-  minicube[x, y, v] = t -- just what vectorify does
  ;-  indcube[x, y, v] = cubeindex
  ;-  minicube[x, y, v] = data[cubeindex]
  vectorify, minicube, x = x, y = y, v = v, t = t, $
             mask = (minicube eq minicube)
  assert, min(minicube[x,y,v] eq t)
  assert, min(indcube[x,y,v] eq cubeindex)
  assert, min(minicube[x,y,v] eq data[cubeindex])

  ;- generate the dendrogram
  generate_dendrogram, merger, clusters = clusters, height = height, xlocation = xlocation, $
                       leafnodes = leafnodes

  ; Patch up roundoff error.
  if ~keyword_set(fast) then begin
     mld=abs(median(levels-shift(levels,-1)))
     matchtol=1e-2*mld
     for jj = 0, n_elements(height)-1 do begin
        diffs = min(abs(height[jj]-levels), ind)
        if diffs lt matchtol then height[jj] = levels[ind]      
     endfor
  endif

  ; This reorders the merger matrix so that it corresponds to the order
  ; that the nodes are plotted in the dendrogram.
  newmerger = merger
  newmerger[*] = 0
  for k = 0, n_elements(leafnodes)-1 do begin
     for j = 0, k do begin
        newmerger[k, j] = merger[leafnodes[k], leafnodes[j]]
     endfor
  endfor
  newmerger = newmerger > transpose(newmerger)
  order = leafnodes
  sz = size(minicube)    
  
  cluster_label = labelclusters(height, clusters, $
                                kernels, levels, minicube, $
                                all_neighbors = all_neighbors, fast = fast, $
                                contour_res = contour_res)
  ;-next line is SLOW, but a useful bug checker
;  validate_label, minicube, kernels, cluster_label, clusters 
  cluster_label = cluster_label[x, y, v] ;- vectorify this cube
  
  ;- create a histogram (with reverse indices) of the cluster labels. 
  h = histogram(cluster_label, min = 0, max = n_elements(height), $
                reverse_indices = ri)

  ; Create a topology structure to contain all the information about the
  ; cloud's analysis
  structure = {merger:merger, cluster_label:cluster_label, $
               cluster_label_h: h, cluster_label_ri: ri, $
               levels:levels, clusters:clusters, height:height, $
               kernels:kernels, order:order, $
               newmerger:newmerger, x:x, y:y, v:v, t:t, sz:sz, $
               cubeindex:cubeindex, szdata:szdata, $
               all_neighbors:all_neighbors, xlocation:xlocation, $
               npix : keyword_set(fast) ? npix : -1, $
               fast : keyword_set(fast)}
  
  if arg_present(pointer) then pointer = ptr_new(structure)
end

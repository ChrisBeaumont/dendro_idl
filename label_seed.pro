;+
; PURPOSE:
;  This procedure creates a mask consisting of all pixels greater than
;  a threshhold, connected to a seed pixel
;
; INPUTS:
;  data: A 2- or 3D array
;  thresh: The intensity threshhold. All pixels in the final have
;          intnsities >= thresh
;  seed: A 2- or 3-element array defining the coordinates of the seed
;             value. It is assumed that the seed pixel intensity is >=
;             thresh. 
; OUTPUTS:
;  result: The mask of the connected region of pixels with intensity >
;          thresh, containing seed
;
; KEYWORD PARAMETERS:
;  EXTERNAL: If set, then this will can an external C program to do
;            all of the work. This is faster if the final mask is a
;            small subset of the original data
;  ALL_NEIGHBORS: Set to consider all 8/26 adjacent pixels/voxels as
;                 "connected" to a given location. Else, only the 4/8
;                 pixels/voxels sharing an edge/face are considered
;                 connected.
;
; PROCEDURE:
;  label_seed will try to call an external C program in the file
;  'libdendro.so' if /EXTERNAL is set. You must compile this yourself,
;  from the libdendro.c source code. The MAKE_DLL program can help
;  with this. If IDL cant find libdendro.so, it skips the C program
;  altogether. 
;
; MODIFICATION HISTORY:
;  Feb 2010: Written by Chris Beaumont
;- 
pro label_seed, data, thresh, seed, result, external = external, $
                all_neighbors = all_neighbors

  ndim = size(data, /n_d)
  sz = size(data)

  if n_elements(result) ne n_elements(data) then result = byte(data*0)
  result[*] = 0B

  ;- the c program requires its input in double. Make it so
  is_dbl = size(data,/type) eq 5
  if ~is_dbl && keyword_set(external) then data = double(data)

  lib = file_search('libdendro.so', count = ct)
  if keyword_set(external) && ct gt 0 then begin
     junk = call_external(lib[0], 'fill', $
                          data, result, long(ndim), $
                          long(sz[1]), long(sz[2]), $
                          ndim eq 3 ? long(sz[3]) : 1L, $
                          long(seed[0]), long(seed[1]), $
                          ndim eq 3 ? long(seed[2]) : 0L, $
                          double(thresh), long(all_neighbors),/unload)
  endif else begin
     r = label_region(data gt thresh, all_n = all_neighbors, /ulong)
     if ndim eq 2 then result = (r eq r[seed[0], seed[1]]) $
     else result = (r eq r[seed[0], seed[1], seed[2]])
  endelse
  return
end

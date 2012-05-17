;+
; PURPOSE:
;  This function locates all of the local maxima in a 2- or 3D data
;  cube.
;
; INPUTS:
;  cube: The data cube
;
; KEYWORD PARAMETERS:
;  indcube: If cubein is some subset of a larger cube, than indcube
;           gives the indices of that larger cube corresponding to
;           each pixel in cubein. If this is provided, then the local maxima
;           will be returned in big-cube-coordinates 
;  friends: Local maxima must be greater than all neighbors within a
;           box of size 2*friends+1 x 2*friends+1 in the first 2
;           dimensions. Defaults to 1
;  specfriends: If cubein is 3D, then local maxima must be greater
;               than pixels within a box of size 2xspecfriends+1 in the
;               3rd dimension. Defaults to 1
;  minval: If provided, then local maxima must be greater than minval.
;
; OUTPUTS:
;  A vector of 1D indices to local maxima. Each maxima satisfies the
;  following crieteria:
;    data[max] > data[cen + (dx, dy, dz)], 
;      -friends <= dx <= friends
;      -friends <= dy <= friends
;      -specfriends <= dz <= specfriends
;    none of the neighbors are NAN
;    data[max] > minval
;
; PROCEDURE:
;  This is an IDL wrapper to a C program, which can identify local
;  maxima much more efficiently than IDL can. You must have created
;  the library 'libdendro.so' from the source code libdendro.c
;  for this to work. The procedure MAKE_DLL can help with this. 
;  If IDL can't find this library, it defaults to the (slower) 
;  alllocmax program
;
; MODIFICATION HISTORY:
;  Feb 2010: Written by Chris Beaumont
;  May 2010: Re-added a line of code to define lmaxcube. cnb.
;  May 24, 2010: Changed how libdendro.so is searched for. cnb.
;  Nov 19, 2010: Added minval keyword. cnb.
;  December 14, 2010: Tweaked libdendro.c. Local maxima can no longer
;                     be near nan values.
;-  
function cnb_alllocmax, cube, indcube = indcube, $
                        friends = friends, $
                        specfriends = specfriends, $
                        minval = minval

 ; GET THE SIZE AND MAKE AN EMPTY CUBE TO FLAG WHERE LOCAL MAXIMA
 ; EXIST.
  sz = size(cube)
  ndim = sz[0]
  if sz[0] eq 2 then begin
     specfriends = 0
     sz[3] = 1
  endif
  lmaxcube = bytarr(sz[1], sz[2], sz[3])
; INITIALIZE THE DEFAULT BOX SIZE TO BE +/- ONE PIXEL (3 x 3 BOX)
  if (n_elements(friends) eq 0) then $
    friends = 1
  if n_elements(specfriends) eq 0 then specfriends = 1 

  ;- default minval
  if n_elements(minval) eq 0 then minval=-!values.d_infinity else $
     minval = double(minval)

  lib = find_libdendro()
  if strlen(lib) eq 0 then begin
     message, /continue, 'Could not find libdendro.so library. Defaulting to idl code'
     message, /continue, 'Try generating libdendro.so with MAKE_DLL'
     result = alllocmax(cube, indcube = indcube, friends = friends, $
                        specfriends = specfriends)
     return, result
  endif
 
  ;- 2d case
  if ndim eq 2 then begin     
     case size(cube, /tname) of 
        'DOUBLE': junk = call_external(lib[0], 'alllocmax_2d_double', $
                                       cube, $
                                       long(friends), $
                                       minval, $
                                       long(sz[1]), $
                                       long(sz[2]), $
                                       lmaxcube)
        'FLOAT':junk = call_external(lib[0], 'alllocmax_2d_float', $
                                     cube, $
                                     long(friends), $
                                     minval, $
                                     long(sz[1]), $
                                     long(sz[2]), $
                                     lmaxcube)
        else: junk = call_external(lib[0], 'alllocmax_2d_float', $
                                   float(cube), $
                                   long(friends), $
                                   minval, $
                                   long(sz[1]), $
                                   long(sz[2]), $
                                   lmaxcube)
     endcase
  endif else begin
                                ;- 3d case
     case size(cube, /tname) of 
        'DOUBLE': junk = call_external(lib[0], 'alllocmax_3d_double', $
                                       cube, $
                                       long(friends), $
                                       long(specfriends), $
                                       minval, $
                                       long(sz[1]), $
                                       long(sz[2]), $
                                       long(sz[3]), $
                                       lmaxcube)
        'FLOAT':junk = call_external(lib[0], 'alllocmax_3d_float', $
                                     cube, $
                                     long(friends), $
                                     long(specfriends), $
                                     minval, $
                                     long(sz[1]), $
                                     long(sz[2]), $
                                     long(sz[3]), $
                                     lmaxcube)
        else: junk = call_external(lib[0], 'alllocmax_3d_float', $
                                   float(cube), $
                                   long(friends), $
                                   long(specfriends), $
                                   minval, $
                                   long(sz[1]), $
                                   long(sz[2]), $
                                   long(sz[3]), $
                                   lmaxcube)
     endcase
  endelse

  lmaxind = where(lmaxcube eq 1B, num)
    if (num eq 0) then begin
    message, 'No true local max found, defaulting to high point in data.', /con
    dummy = max(lmaxcube, lmaxind, /nan)
  endif
    ;message, /con, 'found '+strtrim(num,2)+' local maxima'
; IF THE INDEX CUBE IS SUPPLIED AND THERE ARE LOCAL MAXIMA THEN
; SUBSTITUTE THE INDICES FROM THE CUBE FOR THE ACTUAL INDICES
    if ((n_elements(indcube) gt 0)) then begin
       lmaxind = indcube[lmaxind]
    endif   
    
    return, lmaxind
 end 

;- a small test program, comparing execution speeds
pro test
  sz = 4
  data = randomu(seed, sz, sz, sz)
  data[1,1,1]=10
  data[1,1,2]=!values.f_nan
  print, 'cnb_alllocmax:'
  t0 = systime(/seconds)
  x1 = cnb_alllocmax(data, friends = 1, specfriends = 1)
  t1 = systime(/seconds)
  print, time2string(t1 - t0)
  print, 'alllocmax:'
  x2 = alllocmax(data, friends = 1, specfriends = 1)
  print, time2string(systime(/seconds) - t1)
  a = bytarr(sz, sz, sz)
  b = a
  a[x1] = 1 & b[x2] = 1
  help, where(a ne b)
  print, 'x1'
  print, array_indices(data, x1)
  print, 'x2'
  print, array_indices(data, x2)

  print, data
end

;+
; PURPOSE:
;  This function locates all of the local maxima in a 2- or 3D data
;  cube.
;
; INPUTS:
;  cubein: The data cube
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
;-  
function cnb_alllocmax, cubein, indcube = indcube, $
                        friends = friends, $
                        specfriends = specfriends

 ; GET THE SIZE AND MAKE AN EMPTY CUBE TO FLAG WHERE LOCAL MAXIMA
 ; EXIST.
  sz = size(cubein)
  ndim = sz[0]
  if sz[0] eq 2 then begin
     specfriends = 0
     sz[3] = 1
  endif

; INITIALIZE THE DEFAULT BOX SIZE TO BE +/- ONE PIXEL (3 x 3 BOX)
  if (n_elements(friends) eq 0) then $
    friends = 1
  if n_elements(specfriends) eq 0 then specfriends = 1 

 ; SET THE NONSENSICAL INDICES TO NOT-A-NUMBERS
  badind = where(cubein ne cubein, badct)
  cube = cubein
  if (badct gt 0) then begin
    cube[badind] = min(cubein, /nan)
 endif

  lib = file_search('libdendro.so', count = ct)
  if ct eq 0 then begin
     message, /continue, 'Could not find libdendro.so library. Defaulting to idl code'
     message, /continue, 'Try generating libdendro.so with MAKE_DLL'
     return, alllocmax(cubein, indcube = indcube, friends = friends, $
                       specfriends = specfriends)
  endif
 
  ;- 2d case
  if ndim eq 2 then begin     
     junk = call_external(lib[0], 'alllocmax_2d', $
                          double(cube), $
                          long(friends), $
                          long(sz[1]), $
                          long(sz[2]), $
                          lmaxcube)
  endif else begin
  ;- 3d case
     junk = call_external(lib[0], 'alllocmax_3d', $
                          double(cube), $
                          long(friends), $
                          long(specfriends), $
                          long(sz[1]), $
                          long(sz[2]), $
                          long(sz[3]), $
                          lmaxcube)
  endelse

  lmaxind = where(lmaxcube eq 1B, num)
    if (num eq 0) then begin
    message, 'No true local max found, defaulting to high point in data.', /con
    dummy = max(lmaxcube, lmaxind, /nan)
  endif

; IF THE INDEX CUBE IS SUPPLIED AND THERE ARE LOCAL MAXIMA THEN
; SUBSTITUTE THE INDICES FROM THE CUBE FOR THE ACTUAL INDICES
    if ((n_elements(indcube) gt 0)) then begin
       lmaxind = indcube[lmaxind]
    endif   
    
    return, lmaxind
 end 

;- a small test program, comparing execution speeds
pro test
  sz = 100
  data = randomu(seed, sz, sz, sz)

  t0 = systime(/seconds)
  x1 = cnb_alllocmax(data, friends = 5, specfriends = 5)
  t1 = systime(/seconds)
  print, time2string(t1 - t0)
  x2 = alllocmax(data, friends = 5, specfriends = 5)
  print, time2string(systime(/seconds) - t1)
  return
  a = bytarr(sz, sz, sz)
  b = a
  a[x1] = 1 & b[x2] = 1
  help, where(a ne b)
;  help, x1, x2
end

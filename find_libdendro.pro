;+
; PURPOSE:
;  This function searches for the shared library libdendro in the
;  user's IDL path. It tries multiple suffix options (.so, .exe, etc)
;  to attempt platform-independence. 
;
; RETURNS:
;  The full path to the libdendro shared library, if it exists in the IDL
;  search path. If not, an empty string is returned.
;
; MODIFICATION HISTORY:
;  May 24 2010: Written by Chris Beaumont
;-
function find_libdendro
  result = ''
  suffix = ['.so', '.a', '.sl', '.exe', '.dll']
  nsuff = n_elements(suffix)
  for i = 0, nsuff - 1, 1 do begin
     result=file_which('libdendro'+suffix[i])
     if strlen(result) ne 0 then break
  endfor
  return, result
end


;+
; PURPOSE:
;  This procedure generates a dendrogram from a merger matrix. It makes the
;  appropriate call to cluster_tree and dendrogram_mod, and then patches up
;  the results so that the units match those in the merger matrix.
;
; INPUTS:
;  merger: A square, symmetric merger matrix. The value at merger[i,j] denotes
;          the contour level at which two local maxima in a data cube merge
;          with one another.
;
; KEYWORD PARAMETERS:
;  clusters: On output, contains the n_maxima-1 mergers that define the
;  non-leaf branches of the dendrogram. clusters[*,i] lists the 2 structures
;  which merge to create structure [i+n_maxima] in the dendrogram (structures
;  0-n_maxima-1 are leaves, and correspond to the local maxima (i.e. diagonal
;  entries in merger)
;  xlocation: The (arbitrary) x position of each dendrogram structure, used
;             for plotting purposes
;  height: The height (that is, contour level) corresponding to each
;          dendrogram structure.
;  leafnodes: some variable returned by dendrogram_mod. I'm not sure
;             what it encodes...
;
; MODIFICATION HISTORY:
;  June 2010: Isolated from topologize into its own function. Code written by
;             Erik Rosolowsky.
;-  
pro generate_dendrogram, merger, clusters = clusters, $
                         height = height, $
                         xlocation = xlocation, $
                         leafnodes = leafnodes

  ; Use this chunk of code to make a distance metric that measures how
  ; close two kernels are to each other.  The "distance" is larger for
  ; kernels that merge at lower contour levels.
  sz = size(merger)
  d = max(merger)-merger
  d[indgen(sz[1]), indgen(sz[1])] = 0

  ; Find clusters using the RSI program.
  clusters = cluster_tree(d, linkdistance)
  
  ; Calculate the linkages between all the clusters.  This is a modified
  ; version that returns the X location and height of clusters that are
  ; plotted in the final visualization.
  DENDROGRAM_MOD, clusters, linkdistance, outverts, outconn, $
                  LEAFNODES = leafnodes, xlocation = xlocation, $
                  height = height
  
  
  ; This turns the height output from dendrogram_mod into the actual
  ; brightness temperature associated with each cluster.
  dvals = merger[indgen(sqrt(n_elements(merger))), $
                 indgen(sqrt(n_elements(merger)))] 
  leafy = where(outverts[1, *] eq 0)
  outverts[1, *] = max(dvals)-outverts[1, *]
  outverts[1, leafy] = dvals[leafnodes[outverts[0, leafy]]]
  
  ; This does the same thing for the height vector.
  height = max(merger)-height
  height[findgen(n_elements(leafnodes))] = merger[$
                                           findgen(n_elements(leafnodes)), $
                                           findgen(n_elements(leafnodes))]
end

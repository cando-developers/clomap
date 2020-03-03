(in-package :lomap-spanning)

(defun choose-root (hash-table)
    (block nil
       (maphash (lambda (k v) (return (values k v)))

(defclass back-span ()
  ((distance-to-root :initarg :distance-to-root :accessor distance-to-root)
   (back-vertex :initarg :back-vertex :accessor back-vertex)))

(defun span-from-node (vertex distance visited back-span-info)
  (setf (gethash vertex visited) t)
  (setf distance (1+ distance))
  ;; First visit all the other edges and create a back-span to the current vertex
  (let ((other-vertices (loop for edge in (lomap:edges vertex)
                              for other-vertex = (if (eq (lomap:vertex1 edge) vertex)
                                                     (lomap:vertex2 edge)
                                                     (lomap:vertex1 edge))
                              do (unless (gethash other-vertex visited) ; if we have NOT visited other-vertex then...
                                   (let ((back-span (make-instance 'back-span
                                                                   :distance-to-root distance
                                                                   :back-vertex vertex)))
                                     (setf (gethash other-vertex back-span-info) back-span)))
                              collect other-vertex)))
    ;; then loop over the other vertices and visit them
    (loop for other-vertex in other-vertices
          do (unless (gethash other-vertex visited)
               (span-from-node other-vertex distance visited back-span-info)))))
        
(defun calculate-spanning-tree (graph root-vertex)
  (let ((visited (make-hash-table))
        (back-span-info (make-hash-table)))
    (span-from-node root-vertex 0 visited back-span-info)
    back-span-info))

(defun edge-in-spanning-tree-p (edge back-span-info)
  ;; Check the first vertex of edge
  (let* ((vertex (vertex1 edge))
         (bsi (gethash vertex back-span-info))
         (back-vertex (back-vertex bsi)))
    (when (eq back-vertex (vertex2 edge))
      (return-from edge-in-spanning-tree-p t)))
  (let* ((vertex (vertex2 edge))
         (bsi (gethash vertex back-span-info))
         (back-vertex (back-vertex bsi)))
    (when (eq back-vertex (vertex1 edge))
      (return-from edge-in-spanning-tree-p t)))
  nil)

(defun edges-outside-of-spanning-tree (graph back-span-info)
  (loop for edge in (edges graph)
        when (not (edge-in-spanning-tree-p edge back-span-info))
          collect edge))


(defun generate-backspan-bitvec (vertex back-span-info graph)
  (let ((num-vertices (length (vertices graph)))
        (bitvec (make-array (bitvec-length num-vertices) :element-type 'bit :adjustable nil :initial-element 0)))
    (loop for vertex-index = (index vertex)
          for bsi = (gethash vertex back-span-info)
          do (when bsi
               (let* ((back-vertex (back-vertex bsi))
                      (back-vertex-index (index back-vertex)))
                 (let ((i-index (min vertex-index back-vertex-index))
                       (j-index (max vertex-index back-vertex-index))
                       (bit-index (index-to-bit i-index j-index num-vertices)))
                   (setf (elt bitvec bit-index) 1)
                   (setf vertex back-vertex)))))
    bitvec))
  

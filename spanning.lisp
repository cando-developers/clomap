(in-package :lomap-spanning)

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

   

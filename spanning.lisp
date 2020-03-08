(in-package :lomap)

(defclass back-span ()
  ((distance-to-root :initarg :distance-to-root :accessor distance-to-root)
   (back-vertex :initarg :back-vertex :accessor back-vertex)))

(defparameter *debug-spanning-tree* nil)
(defparameter *largest-distance* nil)

(defun span-from-node (vertex distance back-span-info)
  (when *debug-spanning-tree*
    (format t "On vertex: ~a~%" vertex))
  (setf distance (1+ distance))
  (setf *largest-distance* (max *largest-distance* distance))
  ;; First visit all the other edges and create a back-span to the current vertex
  (let ((other-vertices nil))
    (loop for edge in (edges vertex)
          for other-vertex = (if (eq (vertex1 edge) vertex)
                                 (vertex2 edge)
                                 (vertex1 edge))
          do (multiple-value-bind (value foundp)
                 (gethash other-vertex back-span-info)
               (unless foundp
                 (when *debug-spanning-tree*
                   (format t "   child: ~a~%" other-vertex))
                 (let ((back-span (make-instance 'back-span
                                                 :distance-to-root distance
                                                 :back-vertex vertex)))
                   (setf (gethash other-vertex back-span-info) back-span))
                 (push other-vertex other-vertices))))
    ;; then loop over the other vertices and visit them
    (loop for other-vertex in other-vertices
          do (span-from-node other-vertex distance back-span-info))))

(defun calculate-spanning-tree (graph root-vertex &key debug)
  "Return a hash-table of back-span records that describe a spanning tree from the root.
Also return a second value that is the longest path in the spanning tree from the root."
  (let ((*debug-spanning-tree* debug)
        (*largest-distance* 0))
    (let ((back-span-info (make-hash-table)))
      (setf (gethash root-vertex back-span-info) nil)
      (span-from-node root-vertex 0 back-span-info)
      (values back-span-info *largest-distance*))))

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
  (let* ((num-vertices (length (vertices graph)))
        (bitvec (make-array (bitvec-length num-vertices) :element-type 'bit :adjustable nil :initial-element 0)))
    (loop named walk
          for vertex-index = (index vertex)
          for bsi = (gethash vertex back-span-info)
          do (if bsi
               (let* ((back-vertex (back-vertex bsi))
                      (back-vertex-index (index back-vertex))
                      (bit-index (index-to-bit vertex-index back-vertex-index num-vertices)))
                 (setf (elt bitvec bit-index) 1)
                 (setf vertex back-vertex))
               (return-from generate-backspan-bitvec bitvec)))))

(defun nodes-in-fundamental-cycles-p (graph spanning-tree)
  ;;question: to do check might need to geberate spanning tree here instead of pass one
  ;;becuase this will be a new "chopped graph" when checking constraints.
  (let* ((outedges (edges-outside-of-spanning-tree graph spanning-tree))
         (num-vertices (length (vertices graph)))
         (all-cycles-bitvec (make-array (bitvec-length num-vertices)
                                        :element-type 'bit
                                        :adjustable nil
                                        :initial-element 0)))
    (loop for outedge-i in outedges
          do (let* ((s1 (generate-backspan-bitvec (vertex1 outedge-i) spanning-tree graph))
                    (s2 (generate-backspan-bitvec (vertex2 outedge-i) spanning-tree graph))
                    (s12 (bit-xor s1 s2)))
                    (setf (elt s12 (index-to-bit (index (vertex1 outedge-i)) (index (vertex2 outedge-i)) num-vertices)) 1)
               (setf all-cycles-bitvec (bit-ior s12 all-cycles-bitvec))))
   (format t "~a~%" all-cycles-bitvec)
  (let ((nodes-in-cycles (make-array num-vertices :element-type 't :adjustable nil :initial-element nil)))
    (loop for bit from 0 below (length all-cycles-bitvec)
          do (multiple-value-bind (node1-index node2-index)
                 (bit-to-index bit num-vertices)
               (setf (elt nodes-in-cycles node1-index) t
                     (elt nodes-in-cycles node2-index) t)))
    (format t "~a~%~a~%" nodes-in-cycles (every #'identity nodes-in-cycles)))))

(defun graph-wider-than-p (graph maxdist)
  "Return T if the graph-diameter is less than or equal to maxdist"
  (loop for vertex in (vertices graph)
        do (multiple-value-bind (spanning-tree width)
               (calculate-spanning-tree graph vertex)
             (when (> width maxdist)
               (return-from graph-wider-than-p T))))
  nil)


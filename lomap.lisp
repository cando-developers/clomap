
(in-package :lomap)

(defclass edge ()
  ((vertex1 :initarg :vertex1 :accessor vertex1)
   (vertex2 :initarg :vertex2 :accessor vertex2)
   (sim-score :initarg :sim-score :accessor sim-score)))

(defmethod print-object ((edge edge) stream)
  (print-unreadable-object (edge stream)
    (format stream "LOMAP:EDGE ~a-~a :sim-score ~a"
            (chem:get-name (molecule (vertex1 edge)))
            (chem:get-name (molecule (vertex2 edge)))
            (sim-score edge))))

(defclass vertex ()
  ((molecule :initarg :molecule :accessor molecule)
   (index :initarg :index :accessor index)
   (edges :initarg :edges :initform nil :accessor edges)))

(defmethod print-object ((vertex vertex) stream)
  (print-unreadable-object (vertex stream)
    (format stream "LOMAP:VERTEX molecule: ~a" (chem:get-name (molecule vertex)))))

(defclass graph ()
  ((vertices :initarg :vertices :initform nil :accessor vertices)
   (edges :initarg :edges :initform nil :accessor edges)))

(defun copy-graph-with-edge-removed (graph edge)
  "This returns a copy of the graph with the specified edge removed"
  (let ((new-edges (loop for e in (edges graph)
                         when (not (eq e edge))
                           collect e)))
    (make-instance 'graph :vertices (vertices graph)
                          :edges new-edges)))


;;maximumn common substructure rule
(defun mcsr-similarity-score (mol-a mol-b &key (topological-constraint-theta 2))
  (multiple-value-bind (match delta-a delta-b)
      (molecule-graph.max-clique::compare-molecules mol-a mol-b :topological-constraint-theta topological-constraint-theta)
    (let* ((nmcs (length match))
           (na (+ nmcs (length delta-a)))
           (nb (+ nmcs (length delta-b)))
           (na+nb-2nmcs (- (+ na nb) (* 2 nmcs)))
           (s-mincar (min-common-atom-rule mol-a mol-b))
           (s-ecr (ecr-charge-rule mol-a mol-b))
           (beta -0.1)
           (s-mcsr (exp (* beta na+nb-2nmcs)))
           (s (* s-ecr s-mcsr s-mincar)))
;;         (format t "nmcs -> ~a  na -> ~a  nb -> ~a~%" nmcs na nb)
      s)))

(defun similarity-timing (mol-a molecules &key (topological-constraint-theta 2))
  (loop for mol-b from 0 below (length molecules)
        do (let ((start-time (get-internal-real-time)))
             (mcsr-similarity-score (nth mol-a molecules) (nth mol-b molecules)
                                           :topological-constraint-theta topological-constraint-theta)
             (let ((end-time (get-internal-real-time)))
               (format t "~a~C~a~%" (number-of-heavy-atoms (nth mol-b molecules)) #\tab (/ (float (- end-time start-time)) internal-time-units-per-second))))))
 
(defun ecr-charge-rule (mol-a mol-b)
  (let ((charge-a 0.0))
    (Cando:do-atoms (atm mol-a)
      (incf charge-a (chem:get-charge atm)))
    (let ((charge-b 0.0)
          (ecr 0))
      (cando:do-atoms (atm mol-b)
        (incf charge-b (chem:get-charge atm)))
      (cond ((/= charge-a charge-b) (setq ecr 0))
            (t (setq ecr 1)))
      ecr)))

(defun calc-mol-charge (molecule)
  (let ((charge-a 0.0))
    (Cando:do-atoms (atm molecule)
      (incf charge-a (chem:get-charge atm)))
    charge-a))

(defparameter *ths* 7) 

(defun min-common-atom-rule (mol-a mol-b)
  (let ((num-a (number-of-heavy-atoms mol-a))
        (num-b (number-of-heavy-atoms mol-b))
        (mincar 0))
    (cond ((and (< num-a *ths*) (< num-b *ths*)) (setq mincar 0))
          (t (setq mincar 1)))
  mincar))

(defun similarity-matrix (molecules)
  (let ((matrix (make-array (list (length molecules) (length molecules)) :element-type 't)))
    (loop for yindex from 0 below (1- (length molecules))
          for mol-a = (elt molecules yindex)
          do (loop for xindex from (1+ yindex) below (length molecules)
                   for mol-b = (elt molecules xindex)
                   for similarity = (mcsr-similarity-score mol-a mol-b)
                   do (setf (aref matrix xindex yindex) similarity)
                   do (setf (aref matrix yindex xindex) similarity)))
    matrix))

(defun similarity-graph (molecules matrix)
  (let ((vertices (loop for mol in molecules
                        for index from 0
                        collect (make-instance 'vertex :molecule mol :index index))))
    (let ((graph (make-instance 'graph :vertices vertices)))
      (loop for moly below (1- (length molecules))
            for vertexy = (elt vertices moly)
            do (loop for molx from (1+ moly) below (length molecules)
                     for vertexx = (elt vertices molx)
                     for similarity = (aref matrix molx moly)
                     do (when (> similarity 0.5)
                          (let ((edge (make-instance 'edge :vertex1 vertexx
                                                           :vertex2 vertexy
                                                           :sim-score similarity)))
                            (push edge (edges vertexx))
                            (push edge (edges vertexy))
                            (push edge (edges graph))))))
      graph)))

(defun number-of-heavy-atoms (molecule)
 (let ((count 0))
   (cando:do-atoms (atm molecule)
                   (when (/= (chem:get-atomic-number atm) 1)
                     (incf count)))
   count))

;; functions to convert two zero based indecies (i,j) of a 2d (n x n) upper triangle matrix into a matrix element zero based index (k)
(declaim (inline index-to-bit))
(defun index-to-bit (i j n)
  (declare (fixnum i j n))
  (if (> i j) (rotatef i j)) ;; macro to swap i and j
  (let ((k (- (- (+ (- (* n (/ (- n 1) 2)) (* (- n i) (/ (- (- n i) 1) 2))) j) i) 1)))
  k))

(declaim (inline bit-to-index))
(defun bit-to-index (k n)
  (declare (fixnum i j n))
  (let* ((i (- (- n 2) (truncate (- (/ (sqrt (- (+ (* -8 k) (* n (* 4 (- n 1)))) 7)) 2) 0.5))))
         (j (+ (- (+ k i 1) (/ (* n (- n 1)) 2)) (* (- n i) (/ (- (- n i) 1) 2)))))
         (values i j)))

(defun bitvec-length (n)
  "Return the length of the edge bitvec given the number of nodex"
  (/ (* (1- n) n) 2))



(defun edge-bitvec-from-edges (graph)
  (let* ((num-vertices (length (vertices graph)))
         (bitvec-length (bitvec-length num-vertices))
         (bitvec (make-array bitvec-length
                             :element-type 'bit
                             :adjustable nil
                             :initial-element 0)))
    (format t "bitvec-length -> ~a~%" bitvec-length)
    (loop for edge in (edges graph)
          for i = (index (vertex1 edge))
          for j = (index (vertex2 edge))
          for idx = (index-to-bit i j num-vertices)
;;;          do (format t "i j n .. idx -> ~a ~a ~a .. ~a~%" i j num-vertices idx)
          do (setf (elt bitvec idx) 1)
          )
    bitvec))
        
(defun edges-sorted-by-similarity (graph)
  (let* ((edges-copy (copy-list (edges graph)))
         (sorted-edges (sort edges-copy #'< :key #'sim-score)))
    sorted-edges))

  
(defun lomap-graph (graph &key debug (max-width 8))
  (let ((done nil)
        (num-vertices (length (vertices graph)))
        (sorted-edges (edges-sorted-by-similarity graph))
        )
    (let* ((filename (format nil "/tmp/graph0.dot")))
      (format t "Writing graph to ~a~%" filename)
      (draw-graph-to-file graph filename))
    (loop
      for edge in sorted-edges
      for count from 1
      until done
      do (format t "Trying to remove ~a~%" edge)
      do (let* ((new-graph (copy-graph-with-edge-removed graph edge))
                (_ (loop for edge in new-graph
                         do (format t "new-graph edge: ~a~%" edge)))
                (new-bitvec (edge-bitvec-from-edges new-graph))
                (_ (loop for bit-index below (length new-bitvec)
                         for bit-val = (elt new-bitvec bit-index)
                         do (when (= bit-val 1)
                              (multiple-value-bind (node1-index node2-index)
                                  (bit-to-index bit-index num-vertices)
                                (let ((vertex1 (elt (vertices graph) node1-index))
                                      (vertex2 (elt (vertices graph) node2-index)))
                                  (format t "new-bitvec edge ~a - ~a~%" vertex1 vertex2))))))
                (spanning-tree (calculate-spanning-tree new-graph (first (vertices new-graph))))
                (_ (maphash (lambda (vertex backspan)
                              (format t "spanning-tree ~a -> ~a~%" vertex (back-vertex backspan)))
                            spanning-tree))
                )
           (format t "new-graph vertices: ~a edges: ~a~%" (length (vertices new-graph)) (length (edges new-graph)))
           (progn
             (when debug
               (let* ((filename (format nil "/tmp/graph~a.dot" count)))
                 (format t "Writing graph to ~a~%" filename)
                 (draw-span-graph-to-file new-graph filename :edge-bitvec new-bitvec
                                                             :back-span-info spanning-tree
                                                             :extra-edges (list edge))))
             (let ((too-wide (graph-wider-than-p new-graph max-width))
                   (not-all-in-cycles (not (all-nodes-in-fundamental-cycles-p new-graph spanning-tree))))
               (format t "too-wide -> ~a  not-all-in-cycles -> ~a~%" too-wide not-all-in-cycles)
               (setf done (or too-wide not-all-in-cycles))
               (unless done
                 (setf graph new-graph))))
           ))
    ;; Graph is the simplest graph
    graph))


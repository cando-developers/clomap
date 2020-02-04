
(in-package :lomap)

(defclass edge ()
  ((vertex1 :initarg :vertex1 :accessor vertex1)
   (vertex2 :initarg :vertex2 :accessor vertex2)
   (weight :initarg :weight :accessor weight)))

(defclass vertex ()
  ((molecule :initarg :molecule :accessor molecule)
   (edges :initarg :edges :accessor edges)))

(defclass graph ()
  ((vertices :initarg :vertices :accessor vertices)
   (edges :initarg :edges :accessor edges)))

;;`edges` is a hash (list v1 v2) -> (cons (list v1 v2) weight).
;;`vertices` is a hash v -> (list (cons (list v1 v2) weight) (cons (list v3 v4) weight2) ...)

(defmethod graph-instance :after ((g graph) &key vertices edges &allow-other-keys)
  (with-slots ((v-list vertices) (e-list edges)) g
    (dolist (v vertices)
    ;;set up values for vertecies using v as var for each vertex and assign nil
      (setf (gethash v v-list) nil))
    ;;now build an edge from a -> b
    (dolist (e edges)
    ;;we want `edges` to be a hash e (list v1 v2) -> (cons (list v1 v2) weight)  so each edge (value) points to the verticies its connnected to together with a weight.
      (let* ((a (first (car e))) (b (second (car e))) (weight (or (cdr e) 1)) (value (cons (list a b) weight)))
    ;;now hash through and provide, for each vertex, if it is connected to an edge and what edges it is connected to.
       (multiple-value-bind (v1-edge v1-existsp) (gethash a v-list)
         (multiple-value-bind (v2-edge v2-existsp) (gethash b v-list)
           (setf (gethash a v-list) (cons value v1-edge)
                 (gethash b v-list) (cons value v2-edge)
                 (gethash (car value) e-list) value)))))))

(defmethod are-they-connected ((g graph) v1 v2)
  (or (gethash (list v1 v2) (edges g)) (gethash (list v2 v1) (edges g))))




;;naximumn common substructure rule
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
  (let ((vertices (loop mol in molecules
                        collect (make-instance 'vertex :molecule mol))))
    (let ((graph (make-instance 'graph :vertices vertices)))
      (loop for moly below (1- (length molecules))
            for vertexy = (elt vertices moly)
            do (loop for molx from (1+ moly) below (length molecules)
                     for vertexx = (elt vertices molx)
                     for similarity = (aref matrix molx moly)
                     do (when (> similarity 0.1)
                          (let ((edge (make-instance 'edge :vertex1 vertexx
                                                           :vertex2 vertexy
                                                           :weight similarity)))
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

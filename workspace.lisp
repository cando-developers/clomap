(progn
  (asdf:load-system :lomap)
  (format t "Done~%"))

(in-package :lomap)

(setf *default-pathname-defaults* #P"/users/yonkunas/Development/clomap/")

(asdf:load-asd "/users/yonkunas/Development/clomap/lomap.asd")

(defparameter *mols* (sdf:load-sdf-as-list-of-molecules "/Users/yonkunas/Development/fep-benchmark/eg5/ligands.sdf"))

(defparameter *mols11* (subseq *mols* 6 16))

(defparameter *mat11* (similarity-matrix *mols11*))

(defparameter *graph* (similarity-graph *mols11* *mat11*))

(draw-graph-to-file *graph* "/tmp/tmp.dot")

(defparameter *spanning-tree*
  (calculate-spanning-tree
   *graph*
   (first (vertices *graph*))
   :debug t))


(second (vertices *graph*))
(defparameter *span-bitvec* (generate-backspan-bitvec
                             (second (vertices *graph*))
                             *spanning-tree*
                             *graph*))
(print *span-bitvec*)

(draw-graph-to-file *graph* "/tmp/edge.dot" :edge-bitvec *span-bitvec*)

(draw-span-graph-to-file *graph* "/tmp/edge-span.dot" :edge-bitvec *span-bitvec* :back-span-info *spanning-tree*)

(defparameter *outedges* (edges-outside-of-spanning-tree *graph* *spanning-tree*))

(defparameter *e* (first *outedges*))
(defparameter *s1* (generate-backspan-bitvec
                    (vertex1 *e*)
                    *spanning-tree*
                    *graph*))
(defparameter *s2* (generate-backspan-bitvec
                    (vertex2 *e*)
                    *spanning-tree*
                    *graph*))
(defparameter *s12* (bit-xor *s1* *s2*))
(setf (elt *s12* (index-to-bit (index (vertex1 *e*))
                                      (index (vertex2 *e*))
                                      (length (vertices *graph*)))) 1)
;; Now *s12* has a fundamental cycle
(draw-span-graph-to-file *graph* "/tmp/s12.dot" :edge-bitvec *s12* :back-span-info *spanning-tree*)

(nodes-in-fundamental-cycles-p *graph* *spanning-tree*)
(graph-diameter-p *graph*)

(defparameter *a* (make-array 100 :element-type 'bit))
(setf (elt *a* 32) 1)
(setf (elt *a* 31) 1)
(defparameter *b* (make-array 100 :element-type 'bit))


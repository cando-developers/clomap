
(setf *default-pathname-defaults* #P"/users/yonkunas/Development/clomap/")

(asdf:load-asd "/users/yonkunas/Development/clomap/lomap.asd")

(asdf:load-system :lomap)

(defparameter *mols* (sdf:load-sdf-as-list-of-molecules "/Users/yonkunas/Development/fep-benchmark/eg5/ligands.sdf"))

(defparameter *mols11* (subseq *mols* 6 16))

(defparameter *mat11* (lomap::similarity-matrix *mols11*))

(defparameter *graph* (lomap::similarity-graph *mols11* *mat11*))

(lomap-graphviz::draw-graph-to-file *graph* "/tmp/tmp.dot")

(defparameter *a* (make-array 100 :element-type 'bit))
(setf (elt *a* 32) 1)
(setf (elt *a* 31) 1)
(defparameter *b* (make-array 100 :element-type 'bit))


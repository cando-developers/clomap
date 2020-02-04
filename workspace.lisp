
(setf *default-pathname-defaults* #P"/home/meister/Development/fep-benchmark/eg5/")

(asdf:load-asd "/home/meister/Development/clomap/lomap.asd")

(asdf:load-system :lomap)

(defparameter *mols* (sdf:load-sdf-as-list-of-molecules "ligands.sdf"))

(defparameter *mols8* (subseq *mols* 0 8))

(defparameter *mat8* (lomap::similarity-matrix *mols8*))

(defparameter *graph* (lomap::similarity-graph *mols8* *mat8*))



(setf *default-pathname-defaults* #P"/Users/yonkunas/Development/fep-benchmark/eg5/")

(asdf:load-asd "/Users/yonkunas/Development/lomap/lomap.asd")

(asdf:load-system :lomap)

(defparameter *mols* (sdf:load-sdf-as-list-of-molecules "ligands_original.sdf"))

(defparameter *mols8* (subseq *mols* 0 8))

(defparameter *mat8* (lomap::similarity-matrix *mols8*))

(defparameter *graph* (lomap::similarity-graph *mols8* *mat8*))


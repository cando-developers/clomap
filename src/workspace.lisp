(setf *default-pathname-defaults* (pathname (format nil "~a/Development/clomap/" (ext:getenv "HOME"))))

(asdf:load-asd (pathname (format nil "~a/Development/clomap/lomap.asd" (ext:getenv "HOME"))))

(progn
  (asdf:load-system :lomap)
  (format t "Done~%"))

(in-package :lomap)

(defparameter *mols* (sdf:load-sdf-as-list-of-molecules (format nil "~a/Development/fep-benchmark/eg5/ligands.sdf" (ext:getenv "HOME"))))

;;; Smaller graph
(defparameter *mols4* (subseq *mols* 6 10))
(defparameter *mat4* (similarity-matrix *mols4*))
(defparameter *graph* (similarity-graph *mols4* *mat4*))

;;; Larger graph
(defparameter *mols11* (subseq *mols* 6 16))
(defparameter *mat11* (similarity-matrix *mols11*))
(defparameter *graph* (similarity-graph *mols11* *mat11*))

(graph-set-positions *graph*
                     '(("CHEMBL1096002" 213.29 104.6)
                       ("CHEMBL1084143" 104.4 18)
                       ("CHEMBL1085692" 161.95 47.491)
                       ("CHEMBL1085666" 313.98 194.15)
                       ("CHEMBL1083836" 270.01 156.43)
                       ("CHEMBL1077204" 77.91 73.688)
                       ("CHEMBL1096003" 137.04 100.11)
                       ("CHEMBL1083517" 351.84 163.16)
                       ("CHEMBL1086410" 289.73 110.38)
                       ("CHEMBL1077227" 347.21 114.61)
                       ))
 *graph*

(progn
  (time (defparameter *new-graph* (lomap-graph *graph* :debug nil)))
  (format t "Done~%"))




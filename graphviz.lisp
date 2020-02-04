(in-package :lomap-graphviz)

(defun draw-graph-to-stream (graph stream)
  (format stream "graph goofy {~%")
  (let ((ids (make-hash-table)))
    (loop for vertex in (lomap:vertices graph)
          do (let ((id (gensym)))
               (setf (gethash vertex ids) id)
               (format stream "~a [label = \"~a\"];~%" id (chem:get-name (lomap:molecule vertex)))))
    (loop for edge in (lomap:edges graph)
          do (let ((id1 (gethash (lomap:vertex1 edge) ids))
                   (id2 (gethash (lomap:vertex2 edge) ids)))
               (format stream "~a -- ~a [label = \"~,3f\"];~%" id1 id2 (lomap:weight edge)))))
  (format stream "}~%"))

(defun draw-graph-to-file (graph filename)
  (with-open-file (fout filename :direction :output :if-exists :supersede)
    (draw-graph-to-stream graph fout)))



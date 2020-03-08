(in-package :lomap)

(defun draw-graph-to-stream (graph stream &key edge-bitvec)
  (format stream "graph goofy {~%")
  (let ((ids (make-hash-table))
        (num-vertices (length (vertices graph))))
    (loop for vertex in (vertices graph)
          do (let ((id (gensym)))
               (setf (gethash vertex ids) id)
               (format stream "~a [label = \"~a\"];~%" id (chem:get-name (molecule vertex)))))
    (loop for edge in (edges graph)
          do (let* ((id1 (gethash (vertex1 edge) ids))
                    (id2 (gethash (vertex2 edge) ids))
                    (weight (if edge-bitvec
                                (let* ((index1 (index (vertex1 edge)))
                                       (index2 (index (vertex2 edge)))
                                       (bit-index (index-to-bit index1 index2 num-vertices))
                                       (bit (elt edge-bitvec bit-index)))
                                  (if (= bit 0) 1 3))
                                1)))
               (format stream "~a -- ~a [label = \"~,3f\",penwidth=~a];~%" id1 id2 (sim-score edge) weight))))
  (format stream "}~%"))

(defun draw-span-graph-to-stream (graph stream &key edge-bitvec back-span-info)
  ;;maybe change back-span-info to spanning-tree to make it more readable.
  (format stream "digraph goofy {~%")
  (let ((ids (make-hash-table))
        (num-vertices (length (vertices graph))))
    (loop for vertex in (vertices graph)
          do (let ((id (gensym)))
               (setf (gethash vertex ids) id)
               (format stream "~a [label = \"~a\"];~%" id (chem:get-name (molecule vertex)))))
    (loop for edge in (edges graph)
          do (let* ((id1 (gethash (vertex1 edge) ids))
                    (id2 (gethash (vertex2 edge) ids))
                    (weight (if edge-bitvec
                                (let* ((index1 (index (vertex1 edge)))
                                       (index2 (index (vertex2 edge)))
                                       (bit-index (index-to-bit index1 index2 num-vertices))
                                       (bit (elt edge-bitvec bit-index)))
                                  (if (= bit 0) 1 3))
                                1)))
               ;; Below here id1 and id2 may be swapped
               (let* ((edge-symbol (if back-span-info
                                       (cond
                                         ((let ((back-span (gethash (vertex1 edge) back-span-info)))
                                            (if (and back-span (eq (back-vertex back-span) (vertex2 edge)))
                                                "->"
                                                nil)))
                                         ((let ((back-span (gethash (vertex2 edge) back-span-info)))
                                            (if (and back-span (eq (back-vertex back-span) (vertex1 edge)))
                                                (progn
                                                  (rotatef id1 id2)
                                                  "->")
                                                nil)))
                                         (t nil))
                                       nil))
                      (dir (if edge-symbol
                               "forward"
                               "none"))
                      (edge-symbol (if edge-symbol
                                       edge-symbol
                                       "->")))
                 (format stream "~a ~a ~a [dir=~a,label = \"~,3f\",penwidth=~a];~%"
                         id1 edge-symbol id2 dir (sim-score edge) weight)))))
  (format stream "}~%"))

(defun draw-graph-to-file (graph filename &key edge-bitvec)
  (with-open-file (fout filename :direction :output :if-exists :supersede)
    (draw-graph-to-stream graph fout :edge-bitvec edge-bitvec)))

(defun draw-span-graph-to-file (graph filename &key edge-bitvec back-span-info)
  (with-open-file (fout filename :direction :output :if-exists :supersede)
    (draw-span-graph-to-stream graph fout :edge-bitvec edge-bitvec :back-span-info back-span-info)))

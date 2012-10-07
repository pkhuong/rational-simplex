
;;;; Support

(defun chebyshev-nodes (from to count)
  "Chebyshev nodes tend to result in very close to optimal
   polynomials.  There are solid proofs for smooth functions,
   in fact.

  Enumerate a finite number, but compensate for FP rounding
  by including the immediate neighbourhood of each node's
  rounded value."
  (let ((points (make-array (* 3 (1+ count))
                            :adjustable t :fill-pointer 0))
        (seen   (make-hash-table)))
    (flet ((point (x)
             (let ((x (round-to-float x)))
               (unless (or (gethash x seen)
                           (< x from)
                           (> x to))
                 (setf (gethash x seen) t)
                 (vector-push-extend (make-point x) points)))))
      (point from)
      (point (float-from-bits (float-bits from) 1))
      (point to)
      (point (float-from-bits (float-bits to) -1))
      (loop with stride = (/ (- to from) 2)
            for i from 1 upto count
            for k = (1+ (pull-bits
                         (computable-reals:cos-r (computable-reals:*r
                                                  computable-reals:+pi-r+
                                                  (/ (1- (* 2 i))
                                                     (* 2 count))))))
            for x = (round-to-float (+ from (* k stride)))
            do (point x)
               (point (float-from-bits (float-bits x) -1))
               (point (float-from-bits (float-bits x) 1))))
    (values points seen)))

(defun eval-coefs (points coefs &optional round)
  (declare (optimize debug))
  (let ((poly (make-poly coefs round)))
    (reduce #'max points
            :key (lambda (x)
                   (rational (abs (- (funcall poly (point-loc x))
                                     (point-value x))))))))

(defvar *initial-approximation* (ash 1 8))
(defvar *gap* 5d-2)
(defvar *lp-gap* 1d-4)

(defstruct (node)
  (value sb-ext:double-float-negative-infinity)
  (basis nil)
  (bounds nil))

(defstruct state
  points
  point-set
  primal-solutions
  nodes
  cutoff)


;;;; Root processing
(defun root (from to &key bounds basis points (cutoff 1d300))
  (multiple-value-bind (points seen)
      (if points
          (let ((seen (make-hash-table)))
            (map nil (lambda (x)
                       (setf (gethash x seen) t))
                 points)
            (values (make-array (length points)
                                :adjustable t :fill-pointer (length points)
                                :initial-contents points)
                    seen))
          (chebyshev-nodes from to *initial-approximation*))
    (loop for i upfrom 0 do
     (multiple-value-bind (diff coefs basis)
         (solve-fit points :bounds bounds :advanced-basis basis)
       (multiple-value-bind (actual-diff new-extrema distance root-count)
           (find-error-extrema coefs points)
         (format *trace-output*
                 "Root: ~E ~~ ~E, ~A new extrema (~F bits), ~A roots~%"
                 diff actual-diff (length new-extrema)
                 (if (zerop distance) 0 (log distance 2d0))
                 root-count)
         (when (or (>= (+ diff (* *lp-gap* (abs diff))) actual-diff)
                   (>= actual-diff cutoff)
                   (> i 4))
           (return (values (make-state :points points
                                       :point-set seen
                                       :primal-solutions (make-hash-table :test #'equalp)
                                       :nodes '()
                                       :cutoff cutoff)
                           coefs diff basis)))
         (map nil (lambda (extremum)
                    (let ((extremum (round-to-float extremum)))
                      (unless (gethash extremum seen)
                        (setf (gethash extremum seen) t)
                        (vector-push-extend (make-point extremum)
                                            points))))
              new-extrema))))))


;;;; State update
(defvar *state*)

(defun update-points ()
  (let ((solutions (state-primal-solutions *state*))
        (points    (state-points *state*)))
    (maphash (lambda (k v) v
               (let ((value (eval-coefs points k)))
                 (setf (gethash k solutions) value)))
             solutions))
  (map nil (lambda (node)
             (setf (node-basis node) nil))
       (state-nodes *state*)))

(defun best-upper-bound ()
  (let ((bound sb-ext:double-float-positive-infinity)
        (solution nil))
    (maphash (lambda (k v)
               (when (< v bound)
                 (setf bound v
                       solution k)))
             (state-primal-solutions *state*))
    (values bound solution)))

(defun adjoin-solution (solution &aux (primals (state-primal-solutions *state*))
                                   (points (state-points *state*)))
  (unless (gethash solution primals)
    (setf (gethash solution primals)
          (eval-coefs points solution))))

(defun round-coefs (coefs)
  (map 'simple-vector
       (lambda (x)
         (multiple-value-bind (lo hi)
             (split-float x)
           (if (null lo)
               x
               (let* ((lo (rational lo))
                      (hi (rational hi))
                      (lcm (lcm (denominator lo)
                                (denominator hi)
                                (denominator x))))
                 (if (< (/ (* (- hi lo) (random lcm)) lcm)
                        (- x lo))
                     hi lo)))))
       coefs))


;;; Node exploration
(defun eval-node (bounds basis)
  (if (every (lambda (x)
               (destructuring-bind (lo . hi) x
                 (and lo hi (eql lo hi))))
             bounds)
      (let ((coefs (map 'simple-vector #'car bounds)))
        (values (eval-coefs (state-points *state*) coefs)
                coefs
                nil
                t))
      (solve-fit (state-points *state*) :bounds bounds :advanced-basis basis)))

(defun branch-on-node (coefs bounds value basis)
  (let ((i (position-if #'split-float coefs)))
    (and i
         (let ((coef (elt coefs i))
               (nodes '()))
           (multiple-value-bind (lo hi)
               (split-float coef)
             (format *trace-output*
                     " branch on ~A ~A~%"
                     i (mapcar (lambda (x)
                                 (destructuring-bind (lo . hi) x
                                   (cond ((and lo hi)
                                          (abs (- (float-bits hi)
                                                  (float-bits lo))))
                                         (lo :lo)
                                         (hi :hi)
                                         (t  '-))))
                               bounds))
             (let ((bounds (copy-seq bounds)))
               (setf (elt bounds i) (cons (car (elt bounds i)) lo))
               (push bounds nodes))
             (let ((bounds (copy-seq bounds)))
               (setf (elt bounds i) (cons hi (cdr (elt bounds i))))
               (push bounds nodes))
             (let* ((inf   sb-ext:double-float-positive-infinity)
                    (nodes (sort nodes #'<
                                 :key (lambda (node)
                                        (destructuring-bind (lo . hi)
                                            (elt node i)
                                          (min (abs (- coef (or lo inf)))
                                               (abs (- coef (or hi inf)))))))))
               (mapcar (lambda (node)
                         (make-node :value value
                                    :basis basis
                                    :bounds node))
                       nodes)))))))

(defun node (bounds basis)
  (multiple-value-bind (value coefs basis)
      (eval-node bounds basis)
    (adjoin-solution (round-coefs coefs))
    (cond ((every (lambda (x)
                    (eql x (round-to-float x)))
                  coefs)
           ;; generate cuts
           (multiple-value-bind (actual-diff extrema distance root-count)
               (find-error-extrema coefs (state-points *state*))
             (format *trace-output*
                     "Leaf: ~E ~~ ~E, ~A new extrema (~F bits), ~A roots~%"
                     value actual-diff (length extrema)
                     (if (zerop distance) 0 (log distance 2d0))
                     root-count)
             (cond ((>= (+ value (* *lp-gap* (abs value))) actual-diff)
                    (let ((delta nil))
                      (map nil (lambda (extremum)
                                 (let ((extremum (round-to-float extremum)))
                                   (unless (gethash extremum (state-point-set *state*))
                                     (setf (gethash extremum (state-point-set *state*)) t
                                           delta t)
                                     (vector-push-extend (make-point extremum)
                                                         (state-points *state*)))))
                           extrema)
                      (when delta
                        (update-points))
                      (values actual-diff nil)))
                   (t
                    (values value nil)))))
          (t
           (values value (branch-on-node coefs bounds value basis))))))


;;; Search loop
(defun find-best-node ()
  (let* ((nodes (state-nodes *state*))
         (best (first nodes)))
    (dolist (node (rest nodes) (progn
                                 (setf (state-nodes *state*)
                                       (delete best nodes))
                                 best))
      (when (< (node-value node) (node-value best))
        (setf best node)))))

(defun prune-node (node)
  (let* ((bound (min (best-upper-bound)
                     1d300))
         (bound (min (- bound (* (abs bound) *gap*))
                     (state-cutoff *state*))))
    (and node (< (node-value node) bound) node)))

(defun prune-nodes ()
  (let* ((bound (min (best-upper-bound)
                     (or (state-cutoff *state*)
                         1d300)))
         (bound (min (- bound (* (abs bound) *gap*))
                     (state-cutoff *state*))))
    (setf (state-nodes *state*)
          (delete-if (lambda (node)
                       (>= (node-value node) bound))
                     (state-nodes *state*)))))

(defun branch-and-cut (degree from to function df d2f &key points bounds basis (cutoff 1d300))
  (declare (optimize debug))
  (when bounds
    (assert (= (length bounds) (1+ degree))))
  (let ((*dimension* (1+ degree))
        (*loc-value* function)
        (*loc-dvalue* df)
        (*loc-d2value* d2f)
        *state*
        root-basis
        root-points
        next-node)
    (multiple-value-bind (state coefs diff basis)
        (root from to :bounds bounds :basis basis :points points :cutoff cutoff)
      (setf *state* state
            root-basis basis
            root-points (copy-seq (state-points state)))
      (adjoin-solution (map 'simple-vector #'round-to-float coefs))
      (adjoin-solution (round-coefs coefs))
      (update-points)
      (let ((primal (best-upper-bound)))
        (format *trace-output*
                "At root: ~E ~E (~,3F)~%" diff primal
                (* 100 (/ (- primal diff) primal))))
      (let ((nodes (branch-on-node coefs
                                   (or bounds
                                       (make-list *dimension*
                                                  :initial-element
                                                  '(nil . nil)))
                                   diff basis)))
        (setf next-node (prune-node (pop nodes))
              (state-nodes *state*) nodes)))
    (loop while (or next-node (prune-nodes))
          do (let ((node (or (shiftf next-node nil)
                             (find-best-node))))
               (multiple-value-bind (bound children)
                   (node (node-bounds node) (node-basis node))
                 (when children
                   (setf next-node (prune-node (pop children))
                         (state-nodes *state*) (nconc children (state-nodes *state*))))
                 (prune-nodes)
                 (let* ((primal     (best-upper-bound))
                        (best-bound (reduce #'min (state-nodes *state*)
                                            :key #'node-value
                                            :initial-value (or bound primal))))
                   (format *trace-output* "~A ~E ~E ~E (~,3F)~%"
                           (+ (length (state-nodes *state*)) (if next-node 1 0))
                           (or bound "*") best-bound primal
                           (* 100 (/ (- primal best-bound) primal)))))))
    (multiple-value-call #'values
      (best-upper-bound) *state* root-basis root-points)))

(defun coefs-degree (coefs)
  (1+ (or (position 0 coefs :test-not #'eql :from-end t)
          -1)))

#+nil
(defun enumerate-pareto-front (degree from to function df d2f
                               &optional (cutoff 1d-1)
                                         (cutoff-scale 50d0))
  (declare (optimize debug))
  (let ((primals (make-hash-table :test #'equalp))
        (bound   (make-array (+ degree 2)
                             :initial-element (/ cutoff cutoff-scale)))
        (nogoods '())
        (queue   (make-array 32 :adjustable t :fill-pointer 0))
        (lock    (sb-thread:make-mutex))
        (wait    (sb-thread:make-waitqueue))
        (waitcount 0)
        (done    nil)
        (stem    rational-simplex:*instance-stem*)
        (stdout  *standard-output*))
    (labels ((bound (fixed)
               (branch-and-cut degree from to function df d2f
                               :bounds
                               (loop repeat (1+ degree)
                                     for bound = (pop fixed)
                                     collect (cons bound bound))
                               :cutoff cutoff))
             (rec (fixed skip)
               (format t "solving ~A (~A)~%" fixed (length queue))
               (when (and (= (length fixed) (1+ degree))
                          (not (eql skip :initial)))
                 (let ((fixed fixed))
                   (loop while (and fixed (null (car fixed)))
                         do (pop fixed))
                   (when (every (lambda (x)
                                  (eql x 0))
                                fixed)
                     (format t " pre-eval~%")
                     (return-from rec))))
               (sb-thread:release-mutex lock)
               (when (let ((nfixed (length fixed)))
                       (some (lambda (nogood)
                               (and (<= (length nogood) nfixed)
                                    (every (lambda (nogood fixed)
                                             (or (null nogood)
                                                 (eql nogood fixed)))
                                           nogood fixed)))
                             nogoods))
                 (sb-thread:grab-mutex lock)
                 (format t " nogood~%")
                 (return-from rec))
               (when (eql skip :initial)
                 (setf skip nil))
               (multiple-value-bind (primal coefs state)
                   (if skip
                       (values nil skip nil)
                       (bound fixed))
                 (sb-thread:grab-mutex lock)
                 (unless skip
                   (format t "solved ~A -> ~E" fixed primal)
                   (loop for i from (coefs-degree coefs)
                           below (length bound)
                         do (setf (aref bound i)
                                  (min primal (aref bound i))))
                   (maphash (lambda (k v)
                              (setf (gethash k primals) v))
                            (state-primal-solutions state))
                   (let ((degree (coefs-degree fixed)))
                     (when (plusp degree)
                       (let ((cutoff (* (aref bound (1- degree))
                                        cutoff-scale)))
                         (when (or (>= primal cutoff)
                                   (zerop (hash-table-count
                                           (state-primal-solutions state))))
                           (format t " / pruned~%")
                           (push fixed nogoods)
                           (return-from rec))))
                     (format t "~%")))
                 (let ((fixed fixed))
                   (loop for i from (length fixed) below (length coefs)
                         do (let* ((coef (elt coefs i))
                                   (sign (signum coef)))
                              (loop for j from (min 2 (ceiling (abs coef))) 
                                      downto 0
                                    do
                                       (vector-push-extend
                                        (cons (append fixed (list (* sign j)))
                                              (and (= (* sign j) coef)
                                                   coefs))
                                        queue)))
                            (setf fixed (append fixed '(nil)))))
                 (sb-thread:condition-broadcast wait)))
             (worker ()
               (let ((rational-simplex:*instance-stem*
                       (concatenate 'string stem
                                    (format nil "-~A"
                                            (sb-thread::thread-os-thread
                                             sb-thread:*current-thread*))))
                     (*standard-output* stdout)
                     (*trace-output*    (make-broadcast-stream))
                     computable-reals::(+log2-r+ (+r (ash-r (log-r2 1/7) 1) (log-r2 1/17)))
                     computable-reals::(+PI-R+ (-r (ash-r (atan-r1 1/10) 5)
                                                   (ash-r (atan-r1 1/515) 4)
                                                   (ash-r (atan-r1 1/239) 2)))
                     computable-reals::(+PI/2-R+ (ash-r +pi-r+ -1))
                     computable-reals::(+PI/4-R+ (ash-r +pi-r+ -2)))
                 (sb-thread:with-mutex (lock)
                   (loop until done do
                     (loop while (and (zerop (length queue))
                                      (not done))
                           do (incf waitcount)
                              (sb-thread:condition-wait wait lock)
                              (decf waitcount))
                     (when done (return))
                     (destructuring-bind (fixed . skip)
                         (vector-pop queue)
                       (multiple-value-bind (ok error)
                           (ignore-errors
                            (rec fixed skip)
                            t)
                         (unless ok
                           (unless (eql (sb-thread:mutex-owner lock)
                                        sb-thread:*current-thread*)
                             (sb-thread:grab-mutex lock))
                           (format t "Error: ~A~%" error)))))))))
      (let ((threads (loop repeat 12
                           collect
                           (sb-thread:make-thread #'worker
                                                  :name "worker"))))
        (sb-thread:with-mutex (lock)
          (vector-push-extend (cons '() nil) queue)
          (loop for i from 0 upto degree do
            (vector-push-extend (cons (append
                                       (make-list i)
                                       (make-list (- (1+ degree) i)
                                                  :initial-element 0))
                                      :initial)
                                queue))
          (sb-thread:condition-broadcast wait))
        (unwind-protect
             (progn
               (loop until (sb-thread:with-mutex (lock)
                             (when (and (= waitcount (length threads))
                                        (zerop (length queue)))
                               (setf done t)
                               (sb-thread:condition-broadcast wait)
                               t))
                     do (sleep 60))
               (map nil #'sb-thread:join-thread threads))
          (setf done t)))
      primals)))


(defstruct int-node
  ;; misnomer: it's actually degree+1 (# coefficients...)
  degree bounds value)

(defstruct search-state
  ;; per degree
  nogoods
  pre-exploreds
  bounds
  ;; global
  primals
  queue
  lock
  cvar
  ;; read-only
  points
  cutoff-scale
  from to
  function df d2f)

(defun prune-trailing-nil (x)
  (let ((rev (reverse x)))
    (loop while rev
          do (if (null (car rev))
                 (pop rev)
                 (return)))
    (reverse rev)))

(defun explore-one-node (node state)
  (declare (optimize debug)
           (type search-state state)
           (type int-node node))
  (let* ((degree (int-node-degree node))
         (bounds (int-node-bounds node))
         (value  (int-node-value  node))
         (cutoff (cond ((>= degree 3)
                        (aref (search-state-bounds state)
                              (- degree 3)))
                       ((plusp degree)
                        (aref (search-state-bounds state) 0))
                       (t 1d300)))
         (points (aref (search-state-points state) degree)))
    (assert (= (length bounds) degree))
    (when (and (plusp degree)
               (eql 0 (elt bounds (1- degree))))
      (return-from explore-one-node))
    (unless (some #'null bounds)
      (let* ((key (coerce bounds 'simple-vector))
             (value (eval-coefs points key)))
        (sb-thread:with-mutex ((search-state-lock state))
          (setf (gethash key (search-state-primals state)) value)
          (when (null bounds)
            (let ((cutoffs (search-state-bounds state)))
              (loop for i from degree below (length cutoffs)
                    do (setf (aref cutoffs i)
                             (min (aref cutoffs i) value)))))))
      (return-from explore-one-node))
    (when (and value (> value cutoff))
      (sb-thread:with-mutex ((search-state-lock state))
        (push bounds (aref (search-state-nogoods state) degree)))
      (return-from explore-one-node))
    (let ((nogoods (aref (search-state-nogoods state) degree)))
      (when (some (lambda (nogood)
                    (every (lambda (nogood fixed)
                             (or (null nogood)
                                 (eql nogood fixed)))
                           nogood bounds))
                  nogoods)
        (return-from explore-one-node)))
    (assert (plusp degree))
    (multiple-value-bind (value coefs results)
        (branch-and-cut (1- degree)
                        (search-state-from state)
                        (search-state-to state)
                        (search-state-function state)
                        (search-state-df state)
                        (search-state-d2f state)
                        :bounds
                        (mapcar (lambda (x)
                                  (cons x x))
                                bounds)
                        :cutoff cutoff
                        :points points)
      (sb-thread:with-mutex ((search-state-lock state))
        (when (every #'null bounds)
          (let ((cutoffs (search-state-bounds state)))
            (when value
              (loop for i from degree below (length cutoffs)
                    do (setf (aref cutoffs i)
                             (min (aref cutoffs i) value))))))
        (let ((dst (search-state-primals state)))
          (maphash (lambda (k v)
                     (setf (gethash k dst) v))
                   (state-primal-solutions results))))
      (when (or (null value)
                (zerop (hash-table-count (state-primal-solutions results)))
                (> value cutoff))
        (sb-thread:with-mutex ((search-state-lock state))
          (format t "~A -> ~E / pruned~%" bounds value)
          (let ((bounds (prune-trailing-nil bounds))
                (nogoods (search-state-nogoods state)))
            (loop for i from (max 1 (length bounds)) upto degree do
              (push bounds (aref nogoods i)))))
        (return-from explore-one-node))
      (sb-thread:with-mutex ((search-state-lock state))
        (let ((bounds (prune-trailing-nil bounds))
              #+nil(scale  (search-state-cutoff-scale state))
              (cutoffs (search-state-bounds state))
              (nogoods (search-state-nogoods state)))
          (loop for i from (max 3 (length bounds)) upto degree do
            (if (> value (* #+nil scale (aref cutoffs (- i 3))))
                (push bounds (aref nogoods i))
                (return))))
        (let ((queue (search-state-queue state))
              (seen  (aref (search-state-pre-exploreds state)
                           degree))
              (count 0))
          (flet ((enqueue (i x)
                   (let ((bounds (copy-seq bounds)))
                     (setf (elt bounds i) x)
                     (unless (gethash bounds seen)
                       (setf (gethash bounds seen) t)
                       (incf count)
                       (sb-concurrency:enqueue
                        (make-int-node :degree degree
                                       :bounds bounds
                                       :value  value)
                        queue)))))
            (loop for i upto degree
                  for fixed in bounds
                  for x across coefs do
                    (when (null fixed)
                      (when (and (/= x 0) (< i degree))
                        (enqueue i 0))
                      (let ((abs-x (abs x))
                            (sign  (signum x)))
                        (cond ((eql abs-x 2)
                               (enqueue i sign))
                              ((member abs-x '(0 1 2)))
                              ((< 1 abs-x 3)
                               (enqueue i (* 2 sign))
                               (enqueue i sign))
                              #+nil
                              ((< 1 abs-x 2)
                               (enqueue i (* 2 sign))
                               (enqueue i sign))
                              ((< 0 abs-x 1)
                               (enqueue i sign)))))))
          (sb-thread:condition-broadcast (search-state-cvar state))
          (format t "~A -> ~E (~A / ~A)~%"
                  bounds value (sb-concurrency:queue-count queue) count))))))

(defvar *nthread* 11)

(defun enumerate-pareto-front (degree from to
                               *loc-value* *loc-dvalue* *loc-d2value*
                               &optional (cutoff 1d300) (cutoff-scale 50d0))
  (declare (optimize debug))
  (let ((done    nil)
        (waitcount 0)
        (stem    rational-simplex:*instance-stem*)
        (state   (make-search-state
                  :nogoods (make-array (+ 2 degree)
                                       :initial-element nil)
                  :pre-exploreds (map-into (make-array (+ 2 degree))
                                           (lambda ()
                                             (make-hash-table :test #'equalp)))
                  :bounds  (make-array (+ 2 degree)
                                       :initial-element (/ cutoff cutoff-scale))
                  
                  :primals (make-hash-table :test #'equalp)
                  :queue (sb-concurrency:make-queue)
                  :lock (sb-thread:make-mutex)
                  :cvar (sb-thread:make-waitqueue)
                  
                  :points
                  (coerce (loop for i upto (1+ degree)
                                collect
                                (let ((*dimension* i))
                                  (chebyshev-nodes from to (ash 1 8))))
                          'simple-vector)
                  :from from :to to
                  :cutoff-scale cutoff-scale
                  :function *loc-value*
                  :df *loc-dvalue*
                  :d2f *loc-d2value*))
        (stdout  *standard-output*))
    (flet ((worker (state)
             (let ((rational-simplex:*instance-stem*
                     (concatenate 'string stem
                                  (format nil "-~A"
                                          (sb-thread::thread-os-thread
                                           sb-thread:*current-thread*))))
                   (*standard-output* stdout)
                   (*trace-output*    (make-broadcast-stream))
                   computable-reals::(+log2-r+ (+r (ash-r (log-r2 1/7) 1) (log-r2 1/17)))
                   computable-reals::(+PI-R+ (-r (ash-r (atan-r1 1/10) 5)
                                                 (ash-r (atan-r1 1/515) 4)
                                                 (ash-r (atan-r1 1/239) 2)))
                   computable-reals::(+PI/2-R+ (ash-r +pi-r+ -1))
                   computable-reals::(+PI/4-R+ (ash-r +pi-r+ -2)))
               (loop until done do
                 (sb-thread:with-mutex ((search-state-lock state))
                   (loop while (and (sb-concurrency:queue-empty-p
                                     (search-state-queue state))
                                    (not done))
                         do (incf waitcount)
                            (sb-thread:condition-wait
                             (search-state-cvar state)
                             (search-state-lock state))
                            (decf waitcount)))
                 (when done (return))
                 (loop
                  (multiple-value-bind (node fullp)
                      (sb-concurrency:dequeue (search-state-queue state))
                    (unless fullp
                      (return))
                    (multiple-value-bind (ok error)
                        (ignore-errors
                         (explore-one-node node state)
                         t)
                      (unless ok
                        (sb-thread:with-mutex ((search-state-lock state))
                          (format t "error: ~A~%" error))))))))))
      (let ((threads (loop repeat *nthread*
                           collect
                           (sb-thread:make-thread (lambda ()
                                                    (worker state))
                                                  :name "worker"))))
        (sb-thread:with-mutex ((search-state-lock state))
          (loop for i from 0 upto (1+ degree) do
            (sb-concurrency:enqueue (make-int-node
                                     :degree i
                                     :bounds (make-list i))
                                    (search-state-queue state)))
          (sb-thread:condition-broadcast (search-state-cvar state)))
        (unwind-protect
             (progn
               (loop until (sb-thread:with-mutex ((search-state-lock state))
                             (when (and (= waitcount (length threads))
                                        (sb-concurrency:queue-empty-p
                                         (search-state-queue state)))
                               (setf done t)
                               (sb-thread:condition-broadcast
                                (search-state-cvar state))
                               t))
                     do (sleep 15))
               (map nil #'sb-thread:join-thread threads))
          (setf done t)))
      (values (search-state-primals state)
              state))))

(defun print-hash-table (hash file)
  (with-open-file (s file :direction :output
                          :if-exists :supersede)
    (with-standard-io-syntax
      (write (alexandria:hash-table-alist hash)
             :readably t
             :stream s)))
  t)

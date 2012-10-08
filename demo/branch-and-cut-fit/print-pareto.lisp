(defvar *points*)

(defun eval-poly (coefs)
  (let ((max 0)
        (poly (make-poly coefs t)))
    (map nil (lambda (point)
               (let ((error (abs (- (funcall poly (point-loc point))
                                    (point-value point)))))
                 (setf max (max max error))))
         *points*)
    max))

(defstruct approx
  coefs
  error
  lg-error
  degree
  non-zero
  non-one
  non-two
  constant
  id)

(defun find-non-dominated (approxs)
  (let ((good (make-array (length approxs)
                          :fill-pointer 0)))
    (map nil (lambda (x)
               (unless (or (null x)
                           (some (lambda (y)
                                   (and (not (eql x y)) y
                                        (let ((cmp
                                                (list (- (approx-error y)
                                                         (approx-error x))
                                                      (- (approx-degree y)
                                                         (approx-degree x))
                                                      (- (approx-non-zero y)
                                                         (approx-non-zero x))
                                                      (- (approx-non-one y)
                                                         (approx-non-one x))
                                                      (- (approx-non-two y)
                                                         (approx-non-two x))
                                                      (- (approx-constant y)
                                                         (approx-constant x)))))
                                          (or (and (notany #'plusp cmp)
                                                   (some #'minusp cmp))
                                              (and (every #'zerop cmp)
                                                   (< (approx-id y)
                                                      (approx-id x)))))))
                                 approxs))
                 (vector-push-extend x good)))
         approxs)
    (coerce good 'simple-vector)))

(defun coefs-pareto (coefs from to npoints function)
  (let* ((*loc-value* function)
         (*points* (chebyshev-nodes from to npoints))
         (seen     (make-hash-table :test #'equalp)))
    (find-non-dominated
     (map 'simple-vector
          (lambda (coefs)
            (let* ((coefs (car coefs))
                   (constant (if (plusp (length coefs))
                                 (elt coefs 0)
                                 0))
                   (start (if (plusp (length coefs)) 1 0))
                   (key   (map '(simple-array double-float 1)
                               (lambda (x) (float x 1d0)) coefs)))
              (unless (gethash key seen)
                (setf (gethash key seen) t)
                (let ((error (float (eval-poly coefs) 1d0)))
                  (make-approx
                   :coefs coefs
                   :error error
                   :lg-error (if (zerop error)
                                 sb-ext:double-float-positive-infinity
                                 (floor (- (log error 2d0))))
                   :degree (coefs-degree coefs)
                   :non-zero (count 0 coefs :test-not #'eql
                                    :start start)
                   :non-one  (count-if-not (lambda (x)
                                             (member x '(-1 0 1)))
                                           coefs
                                           :start start)
                   :non-two  (count-if-not (lambda (x)
                                             (member x '(-2 -1 0 1 2)))
                                           coefs
                                           :start start)
                   :constant (case (abs constant)
                               (0 0)
                               (1 1)
                               (2 2)
                               (t 3))
                   :id (hash-table-count seen))))))
          coefs))))

(defun sort-by (values &rest accessors)
  (let ((values (if (simple-vector-p values)
                    (copy-seq values)
                    (coerce values 'simple-vector))))
    (dolist (accessor (reverse accessors) values)
      (setf values (stable-sort values #'< :key accessor)))))

(defvar *accessors* (list (cons #'approx-error "error")
                          (cons (lambda (x)
                                  (- (approx-lg-error x)))
                                "lb_error")
                          (cons (lambda (x)
                                  (1- (approx-degree x)))
                                "degree")
                          (cons #'approx-non-zero "non_zero")
                          (cons #'approx-non-one  "non_one")
                          (cons #'approx-non-two  "non_two")
                          (cons #'approx-constant "constant")))

(defvar *default-accessors* (list (cons (lambda (x)
                                          (- (approx-lg-error x)))
                                        "lb_error")
                                  (cons (lambda (x)
                                          (1- (approx-degree x)))
                                        "degree")
                                  (cons #'approx-error "error")
                                  (cons #'approx-non-zero "non_zero")
                                  (cons #'approx-non-one  "non_one")
                                  (cons #'approx-non-two  "non_two")
                                  (cons #'approx-constant "constant")))

(defvar *default-accessors2* (list (cons (lambda (x)
                                           (1- (approx-degree x)))
                                         "degree")
                                   (cons (lambda (x)
                                           (- (approx-lg-error x)))
                                         "lb_error")
                                   (cons #'approx-non-zero "non_zero")
                                   (cons #'approx-non-one  "non_one")
                                   (cons #'approx-non-two  "non_two")
                                   (cons #'approx-constant "constant")
                                   (cons #'approx-error "error")))


(defun md5 (x)
  #+sbcl
  (let ((sum (sb-md5:md5sum-string x)))
    (format nil "~{~2,'0X~}" (coerce sum 'list)))
  #-sbcl
  (sxhash x))

(defun md5-coefs (coefs)
  #+sbcl
  (ecase *float-mode*
    (double-float
     (let* ((coefs (map '(simple-array double-float 1) (lambda (x)
                                                         (float x 1d0))
                        coefs))
            (nbytes (* 8 (length coefs)))
            (bytes (make-array nbytes :element-type '(unsigned-byte 8))))
       (sb-impl::ub8-bash-copy coefs 0 bytes 0 nbytes)
       (let ((sum (sb-md5:md5sum-sequence bytes)))
         (format nil "~{~2,'0X~}" (coerce sum 'list)))))
    (single-float
     (let* ((coefs (map '(simple-array single-float 1) (lambda (x)
                                                         (float x 1s0))
                        coefs))
            (nbytes (* 4 (length coefs)))
            (bytes (make-array nbytes :element-type '(unsigned-byte 8))))
       (sb-impl::ub8-bash-copy coefs 0 bytes 0 nbytes)
       (let ((sum (sb-md5:md5sum-sequence bytes)))
         (format nil "~{~2,'0X~}" (coerce sum 'list))))))
  #-sbcl
  (sxhash coefs))

(defun print-index (function long-description approx accessors)
  (let ((approx (apply 'sort-by approx (mapcar #'car accessors)))
        (accessors (mapcar #'car accessors))
        (names  (mapcar #'cdr accessors))
        (*read-default-float-format* *float-mode*)
        (error-width  (ecase *float-mode*
                        (single-float 14)
                        (double-float 24))))
    (with-open-file (s (format nil "~A~{-~A~}"
                               function names)
                       :direction :output
                       :if-exists :supersede)
      (format s "~A ~A~%" function long-description)
      (dolist (name names)
        (if (equal name "error")
            (format s "~V,A " error-width "error")
            (format s "~8,A " name)))
      (format s "| coefficients | rationals | hash~%")
      (map nil (lambda (approx)
                 (loop for accessor in accessors
                       for name in names
                       do (let ((x (funcall accessor approx)))
                            (when (equal name "lb_error")
                              (setf x (- x)))
                            (cond ((eql x sb-ext:double-float-positive-infinity)
                                   (format s "~V,A " error-width "inf"))
                                  ((floatp x)
                                   (format s "~V,E " error-width (coerce x *float-mode*)))
                                  (t
                                   (format s "~8,A " x)))))
                 (let ((coefs (approx-coefs approx)))
                   (setf coefs (subseq coefs 0 (coefs-degree coefs)))
                   (format s "| ~{~W ~}| ~{~A ~}| ~A-~A~%"
                           (map 'list (lambda (x)
                                        (floatify x))
                                coefs)
                           (coerce coefs 'list)
                           function
                           (md5-coefs coefs))))
           approx))))

(defun print-all-indices (approx function long-description
                          &key (max-error sb-ext:double-float-positive-infinity))
  (let ((approx (remove-if (lambda (x)
                             (> (approx-error x) max-error))
                           approx))
        (seen (make-hash-table :test #'equalp)))
    (let ((error (car *accessors*)))
      (alexandria:map-permutations
       (lambda (accessors)
         (let* ((split (1+ (position error accessors)))
                (first (subseq accessors 0 split))
                (second (subseq accessors split)))
           (setf second
                 (sort second #'<
                       :key (lambda (x)
                              (position x *accessors*))))
           (let* ((accessors (append first second))
                  (key       (mapcar #'cdr accessors)))
             (unless (gethash key seen)
               (setf (gethash key seen) t)
               (print-index
                function long-description
                approx accessors)))))
       *accessors*))))

(defun print-default-index (approx function long-description
                            &key (max-error
                                  sb-ext:double-float-positive-infinity))
  (let ((approx (remove-if (lambda (x)
                             (> (approx-error x) max-error))
                           approx)))
    (print-index function long-description
                 approx *default-accessors*)
    (print-index function long-description
                 approx *default-accessors2*)))

(defun print-indices ()
  (mapc (lambda (x)
          (destructuring-bind (file from to approximatee
                               name long-description)
              x
            (print-default-index
             (coefs-pareto (with-open-file (s file)
                             (read s))
                           from to (ash 1 13)
                           (coerce approximatee 'function))
             name long-description :max-error 0.1d0)
            (format t "~A done~%" name)
            (force-output)))
        `(
          ("arctan-16" 0 1 ,#'computable-reals:atan-r "atan" "over [0, 1], degree at most 16, error at most 0.1")
          ("arctan-wide-16"  -1 1 ,#'computable-reals:atan-r "wider-atan" "over [-1, 1], degree at most 16, error at most 0.1")
          ("cos-narrow-16" ,(- (float-from-bits (float-bits (/ pi 2)) 2))
                           ,(float-from-bits (float-bits (/ pi 2)) 2)
                           ,#'computable-reals:cos-r
                           "cos" "over [-pi/2, pi/2], degree at most 16, error at most 0.1")
          ("sin-narrow-16" ,(- (float-from-bits (float-bits (/ pi 2)) 2))
                           ,(float-from-bits (float-bits (/ pi 2)) 2)
                           ,#'computable-reals:sin-r
                           "sin" "over [-pi/2, pi/2], degree at most 16, error at most 0.1")
          ("exp-16" 0 1 ,#'computable-reals:exp-r "exp" "over [0, 1], degree at most 16, error at most 0.1")
          ("exp-wide-16" -1 1 ,#'computable-reals:exp-r "wider-exp" "over [-1, 1], degree at most 16, error at most 0.1")
          ("lg1px-16" 0 1 ,(lambda (x)
                             (computable-reals:log-r (computable-reals:+r 1 x) 2))
                      "lg1px" "(log_2 1+x) over [0, 1], degree at most 16, error at most 0.1")
          ("log1px-16" 0 1 ,(lambda (x)
                              (computable-reals:log-r (computable-reals:+r 1 x)))
                       "log1px" "(log 1+x) over [0, 1], degree at most 16, error at most 0.1")
          ("log-16" 1 2 ,(lambda (x)
                           (computable-reals:log-r x))
                    "log" "over [1, 2], degree at most 16, error at most 0.1")
          )))

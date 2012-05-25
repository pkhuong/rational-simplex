(declaim (type (integer 1) *dimension*))
(defvar *dimension* 6)

(defun powers (x &optional (degree (1- *dimension*)))
  (let ((x (rational x)))
    (coerce (loop for i upto degree
                  collect (round-to-double (expt x i)))
            'simple-vector)))

(defvar *loc-parameters* 'powers)
(defvar *loc-value* (lambda (x)
                      (rational
                       (exp (float x 1d0)))))

(defstruct (point
            (:constructor %make-point (loc parameters value)))
  loc
  parameters
  value)

(defun make-point (loc)
  (let ((loc (rational loc)))
    (%make-point loc
                 (funcall *loc-parameters*
                          loc)
                 (funcall *loc-value*
                          loc))))

(defun linexpr-dot (vars values)
  (let ((acc 0))
    (map nil (lambda (var value)
               (lp:addf acc var value))
         vars values)
    acc))

(defun solve-fit (points)
  (lp:with-model (:name "linf fit" :sense :minimize)
    (let ((diff  (lp:var :name "diff" :obj 1))
          (coefs (loop for i below *dimension* collect
                       (lp:var :name (format nil "coef ~A" i)
                               :lower nil))))
      (map nil (lambda (point)
                 (let ((lhs (linexpr-dot coefs
                                         (point-parameters point)))
                       (rhs (point-value point)))
                   (lp:post<= lhs (lp:+ rhs diff))
                   (lp:post>= lhs (lp:- rhs diff))))
           points)
      (multiple-value-bind (status obj values)
          (lp:solve)
        (assert (eql status :optimal))
        (values obj (map 'simple-vector
                         (lambda (var)
                           (gethash var values 0))
                         coefs))))))

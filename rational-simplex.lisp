(defpackage "RATIONAL-SIMPLEX"
  (:use)
  (:nicknames "LP")
  (:export "MODEL" "*MODEL*" "NAME" "SENSE" "OBJ" "CONSTRAINTS"
           "WITH-MODEL"
           "VAR" "LOWER-BOUND" "UPPER-BOUND"
           "LINEAR-EXPRESSION" "LINEXPR" "COEFS" "CONSTANT"
           "CONSTRAINT" "LHS" "CMP" "RHS"

           "ADD" "ADDF" "SUB" "SUBF" "SCALE" "REMOVE-VAR"
           "ADD-CONSTRAINT" "DEL-CONSTRAINT" "CONSTRAIN" "POST"

           "+" "-" "*" "/" "<=" ">=" "="
           "POST<=" "POST>=" "POST="
           
           "%PRINT-MODEL" "PRINT-FLOAT-MODEL" "PRINT-RATIONAL-MODEL"

           "SOLVE" "SOLVE-WITH-PRINTER"
           "*SOLVERS-PATH*" "*INSTANCE-STEM*"))

(defpackage "RATIONAL-SIMPLEX.IMPL"
  (:use "CL" "CL-PPCRE")
  (:import-from "RATIONAL-SIMPLEX"
                "MODEL" "*MODEL*" "NAME" "SENSE" "OBJ" "CONSTRAINTS"
                "WITH-MODEL"
                "VAR" "LOWER-BOUND" "UPPER-BOUND"
                "LINEAR-EXPRESSION" "LINEXPR" "COEFS" "CONSTANT"
                "CONSTRAINT" "LHS" "CMP" "RHS"

                "ADD" "ADDF" "SUB" "SUBF" "SCALE" "REMOVE-VAR"
                "ADD-CONSTRAINT" "DEL-CONSTRAINT" "CONSTRAIN" "POST"
                "POST<=" "POST>=" "POST="
                
                "%PRINT-MODEL" "PRINT-FLOAT-MODEL" "PRINT-RATIONAL-MODEL"

                "SOLVE" "SOLVE-WITH-PRINTER"
                "*SOLVERS-PATH*" "*INSTANCE-STEM*"))

(in-package "RATIONAL-SIMPLEX.IMPL")
(defun maybe-rational (x)
  (if (numberp x)
      (rational x)
      x))
(defun copy-hash-table (x)
  (let ((dst (make-hash-table :test (hash-table-test x))))
    (maphash (lambda (k v)
               (setf (gethash k dst) v))
             x)
    dst))

(defclass model ()
  ((name :reader name  :initarg :name :initform "Unnamed model")
   (sense :reader sense :initarg :sense :initform :minimize)
   (obj  :accessor obj :initarg :obj)
   (constraints :accessor %constraints
                :initarg  :constraints
                :initform (make-hash-table))
   (vars        :reader   %vars
                :initarg  :vars
                :initform (make-hash-table))
   (id-var      :reader id-var
                :initform (make-array 32 :adjustable t
                                         :fill-pointer 0))))

(defgeneric (setf sense) (sense model)
  (:method (sense (model model))
    (check-type sense (member :maximize :minimize))
    (setf (slot-value model 'sense) sense)))

(defgeneric constraints (model))
(defgeneric (setf constraints) (value model))
(defmethod constraints ((model model))
  (copy-hash-table (%constraints model)))
(defmethod (setf constraints) (value (model model))
  (check-type value hash-table)
  (setf (%constraints model) (copy-hash-table value)))

(declaim (type model *model*))
(defvar *model*)

(defclass var ()
  ((name  :reader name  :initarg :name)
   (id    :reader id    :initform (length (id-var *model*)))
   (model :reader %model :initform *model*)
   (lower-bound :reader   lower-bound
                :initarg  :lower-bound
                :initform nil)
   (upper-bound :reader   upper-bound
                :initarg  :upper-bound
                :initform nil)))

(defclass linear-expression ()
  ((coefs :reader   %coefs
          :initarg  :coefs
          :initform (make-hash-table))
   (constant :reader   constant
             :initarg :constant
             :initform 0)))

(declaim (type unsigned-byte *constraint-id*))
(defvar *constraint-id* 0)

(defclass constraint ()
  ((lhs :reader %lhs :initarg :lhs)
   (cmp :reader cmp :initarg :cmp)
   (rhs :reader rhs :initarg :rhs)
   (model :reader %model :initform *model*)
   (id  :reader id :initform (incf *constraint-id*))))

(defun adjust-linearexpr (dst delta-coefs)
  (loop with coefs = (%coefs dst)
        for (var . scale) in delta-coefs do
          (check-type var var)
          (check-type scale real)
          (let ((coef (incf (gethash var coefs 0)
                            (rational scale))))
            (when (zerop coef)
              (remhash var coefs)))
        finally (return dst)))

(defun make-linexpr (&key constant coefs)
  (let ((dst (make-instance 'linear-expression
                            :constant (rational (or constant 0)))))
    (adjust-linearexpr dst coefs)))

(defun copy-linexpr (src &key
                           constant
                           delta-constant
                           delta-coefs)
  (let ((new (make-instance
              'linear-expression
              :coefs (copy-hash-table
                      (%coefs src))
              :constant (rational (+ (or constant
                                         (constant src))
                                     (or delta-constant
                                         0))))))
    (adjust-linearexpr new delta-coefs)
    new))

(defun linexpr (&optional constant &rest coef-var)
  (make-linexpr
   :constant (or constant 0)
   :coefs (loop for (coef var) on coef-var by #'cddr
                do (when (typep coef 'var)
                     (rotatef coef var))
                   (check-type coef real)
                   (check-type var  var)
                collect (cons var coef))))

(defgeneric %add (x y z))
(defun add (x y &optional z)
  (setf z (or z 1))
  (when (realp y)
    (rotatef y z))
  (assert (realp z))
  (%add (maybe-rational x)
        (maybe-rational y)
        (maybe-rational z)))

(defun sub (x y &optional z)
  (setf z (or z 1))
  (add x y (- z)))

(defun scale (x scale)
  (add 0 x scale))

(defun remove-var (linexpr var)
  (check-type linexpr linear-expression)
  (check-type var var)
  (let ((coef (gethash var (%coefs linexpr))))
    (if coef
        (sub linexpr var)
        linexpr)))

(define-modify-macro addf (y &optional (z 1)) add)
(define-modify-macro subf (y &optional (z 1)) sub)
(define-modify-macro scalef (scale) scale)
(define-modify-macro remove-varf (var) remove-var)

(defun rational-simplex:+ (x &rest ys)
  (dolist (y ys x)
    (addf x y)))

(defun rational-simplex:- (x &rest ys)
  (dolist (y ys x)
    (subf x y)))

(defun rational-simplex:* (x y)
  (when (realp x)
    (rotatef x y))
  (check-type y real)
  (scale x y))

(defun rational-simplex:/ (x y)
  (check-type y real)
  (scale x (/ y)))

(defmethod %add ((x real) (y real) z)
  (+ x (* y z)))

(defmethod %add ((x real) (y var) z)
  (make-linexpr :constant x
                :coefs `((,y . ,z))))

(defmethod %add ((x var) (y real) z)
  (make-linexpr :constant (* y z)
                :coefs `((,x . 1))))

(defmethod %add ((x var) (y  var) z)
  (make-linexpr :coefs `((,x . 1)
                         (,y . ,z))))

(defmethod %add ((x linear-expression) (y real)
                 z)
  (copy-linexpr x :delta-constant (* y z)))

(defun linexpr-coefs (linexpr &optional scale)
  (setf scale (or scale 1))
  (let ((coefs '()))
    (maphash (lambda (k v)
               (push (cons k (* v scale))
                     coefs))
             (%coefs linexpr))
    (sort coefs #'< :key (lambda (x)
                           (id (car x))))))

(defmethod %add ((x real) (y linear-expression)
                 z)
  (make-linexpr :constant (+ x (* (constant y) z))
                :coefs (linexpr-coefs y z)))

(defmethod %add ((x linear-expression) (y var)
                 z)
  (copy-linexpr x :delta-coefs `((,y . ,z))))

(defmethod %add ((x var) (y linear-expression)
                 z)
  (make-linexpr :constant (* z (constant y))
                :coefs (acons x 1
                              (linexpr-coefs y z))))

(defmethod %add ((x linear-expression) (y linear-expression)
                 z)
  (copy-linexpr
   x
   :delta-constant (* z (constant y))
   :delta-coefs    (linexpr-coefs y z)))

(defun model (&key name sense)
  (check-type name  (or null string))
  (check-type sense (member nil :minimize :maximize))
  (make-instance 'model :name  (or name "Unnamed model")
                        :obj   (make-linexpr)
                        :sense (or sense :minimize)))

(defmacro with-model ((&key name sense) &body body)
  `(let ((*model* (model :name ,name :sense ,sense))
         (*constraint-id* 0))
     ,@body))

(defun var (&key name (lower 0) upper obj)
  (check-type name  (or null string))
  (check-type lower (or null real))
  (check-type upper (or null real))
  (check-type obj   (or null real))
  (let ((var (make-instance 'var
                            :name name
                            :lower-bound (maybe-rational
                                          lower)
                            :upper-bound (maybe-rational
                                          upper))))
    (when obj
      (addf (obj *model*) var (maybe-rational obj)))
    (setf (gethash var (%vars *model*)) t)
    (let ((id     (id var))
          (id-var (id-var *model*)))
      (assert (= id (length id-var)))
      (vector-push-extend var (id-var *model*)))
    var))

(defmethod print-object ((var var) stream)
  (print-unreadable-object (var stream :type t)
    (format stream "~A [~F:~F]"
            (or (name var)
                (format nil "v~A" (id var)))
            (let ((lower (lower-bound var)))
              (if lower (float lower 1d0) "-inf"))
            (let ((upper (upper-bound var)))
              (if upper (float upper 1d0) "+inf")))))

(defun constraint (lhs cmp rhs)
  (check-type lhs (or var linear-expression real))
  (check-type rhs (or var linear-expression real))
  (check-type cmp (member <= >= =))
  (setf lhs (maybe-rational lhs)
        rhs (maybe-rational rhs))
  (when (and (realp lhs)
             (realp rhs))
    (return-from constraint
      (make-instance 'constraint
                     :lhs '()
                     :cmp cmp
                     :rhs (- rhs lhs))))
  (let* ((linexpr (add lhs rhs -1))
         (coefs   (linexpr-coefs linexpr)))
    (every (lambda (coef)
             (eql *model* (%model (car coef))))
           coefs)
    (make-instance 'constraint
                   :lhs coefs
                   :cmp cmp
                   :rhs (- (constant linexpr)))))

(defgeneric lhs (object))
(defmethod lhs ((constraint constraint))
  (make-linexpr :coefs (%lhs constraint)))

(defun add-constraint (model constraint)
  (check-type constraint constraint)
  (check-type model model)
  (assert (eql model (%model constraint)))
  (setf (gethash constraint (%constraints model)) t)
  constraint)

(defun del-constraint (model constraint)
  (check-type constraint constraint)
  (check-type model model)
  (remhash constraint (%constraints model)))

(defun constrain (lhs cmp rhs &optional model)
  (add-constraint (or model *model*)
                  (constraint lhs cmp rhs)))

(defun rational-simplex:<= (lhs rhs)
  (constraint lhs '<= rhs))
(defun rational-simplex:>= (lhs rhs)
  (constraint lhs '>= rhs))
(defun rational-simplex:= (lhs rhs)
  (constraint lhs '= rhs))

(defun post (lhs cmp rhs &optional model)
  (constrain lhs cmp rhs model))

(defun post>= (lhs rhs &optional model)
  (post lhs '>= rhs model))
(defun post<= (lhs rhs &optional model)
  (post lhs '<= rhs model))
(defun post= (lhs rhs &optional model)
  (post lhs '= rhs model))

(defun print-coefs (coefs numberify stream)
  (format stream "~{~A v~A~^ + ~}"
          (loop for (var . scale) in coefs
                collect (funcall numberify scale)
                collect (id var))))

(defun print-constraints (constraints numberify stream)
  (let ((rows '()))
    (maphash (lambda (k v) v
               (push k rows))
             constraints)
    
    (setf rows (sort rows #'< :key #'id))
    (dolist (row rows)
      (format stream " C~A: " (id row))
      (print-coefs (%lhs row) numberify stream)
      (format stream " ~A ~A~%"
              (ecase (cmp row)
                (<= "<=")
                (>= ">=")
                (=  "="))
              (funcall numberify (rhs row))))))

(defun print-bounds (vars numberify stream)
  (let ((columns '()))
    (maphash (lambda (k v) v
               (push k columns))
             vars)
    (setf columns (sort columns #'< :key #'id))
    (dolist (var columns)
      (let ((lb (lower-bound var))
            (ub (upper-bound var))
            (id (id var)))
        (cond ((or lb ub)
               (format stream " ")
               (when lb
                 (format stream "~A <= "
                         (funcall numberify lb)))
               (format stream "v~A" id)
               (when ub
                 (format stream " <= ~A"
                         (funcall numberify ub)))
               (terpri stream))
              (t
               (format stream " v~A free~%"
                       id)))))))

(defun %print-model (model numberify stream
                     &aux (*read-default-float-format* 'double-float))
  (format stream "Problem~%")
  (format stream " ~A~%" (or (name model)
                             "Unnamed model"))
  (format stream "~A~%"  (ecase (sense model)
                           (:minimize "Minimize")
                           (:maximize "Maximize")))
  (format stream " obj: obj_var~%")
  (format stream "Subject To~%")
  (format stream " obj_con: ")
  (let* ((obj      (obj model))
         (coefs    (linexpr-coefs   obj))
         (constant (constant obj)))
    (when coefs
      (print-coefs coefs numberify stream))
    (format stream " - obj_var = ~A~%"
            (funcall numberify constant)))
  (print-constraints (constraints model) numberify stream)
  (format stream "Bounds~%")
  (print-bounds (%vars *model*) numberify stream)
  (format stream "End~%"))

(defun print-float-model (model stream)
  (let ((*read-default-float-format* 'double-float))
    (%print-model model (lambda (x)
                          (float x 1d0))
                  stream)))

(defun print-rational-model (model stream)
  (let ((*read-default-float-format* 'double-float))
    (%print-model model (lambda (x)
                          (rational x))
                  stream)))

(defvar *solvers-path* (format nil "~A/bin/"
                               (directory-namestring *load-pathname*)))
(defvar *instance-stem* "/tmp/rational-simplex-instance")
(defvar *status-codes* '((1 . :optimal)
                         (2 . :infeasible)
                         (3 . :unbounded)
                         (4 . :unsolved)))

(defun parse-output (s)
  (let ((status      nil)
        (obj         0)
        (var-values  (make-hash-table :test #'equal))
        (time        nil)
        (*read-default-float-format* 'double-float))
    (labels ((parse-status (line)
               (register-groups-bind ((#'parse-integer value))
                   ("LP Value: [-+0-9.]+, status ([0-9])" line)
                 (assert (not status))
                 (setf status value))
               (register-groups-bind ((#'read-from-string sec))
                   ("Time for SOLVER: (.*) seconds." line)
                 (assert (not time))
                 (setf time sec)))
             (parse-value (line)
               (register-groups-bind (name (#'read-from-string value))
                   ("(.*) = ([-+0-9.e/]*)" line)
                 (assert (not (gethash name var-values)))
                 (setf (gethash name var-values) value)))
             (parse-values (stream)
               (loop for line = (read-line stream nil)
                     while line
                     do (parse-value line))))
      (with-input-from-string (s s)
        (loop for line = (read-line s nil)
              while line
              do (parse-status line)
              when (equal line "Solution Values")
                do (return (parse-values s)))))
    (assert (and status time))
    (values obj
            (or (cdr (assoc status *status-codes*))
                (error "Unknown status code ~A" status))
            var-values
            time)))

(defun %solve-with-solver (solver printer basis-p
                           &aux
                             (lp-file (format nil "~A.lp"
                                       *instance-stem*))
                             (bas-file (format nil "~A.bas"
                                        *instance-stem*)))
  (when printer
    (with-open-file (s lp-file
                       :direction :output
                       :if-exists :supersede)
      (funcall printer s)))
  (let ((s #+sbcl (with-output-to-string (s)
                    (sb-ext:run-program (format nil "~A/~A"
                                                *solvers-path* solver)
                                        `("-O"
                                          "-b" ,bas-file
                                          ,@(and basis-p
                                                 `("-B" ,bas-file))
                                          "-L" ,lp-file)
                                        :output s
                                        :error nil
                                        :wait t))
           #-sbcl
           (trivial-shell:shell-command
            (format nil "~A/~A -O -b ~A ~A -L ~A"
                    *solvers-path* solver
                    bas-file (if basis-p
                                 (format nil "-B ~A" bas-file)
                                 "")
                    lp-file))))
    (parse-output s)))

(defun double (x)
  (float x 1d0))

(defvar *solvers* '(("dbl_solver" double)
                    ("ldbl_solver" double)
                    ("float128_solver" double)
                    ("mpf_solver" double)
                    ("mpq_solver" rational)))

(defun %solve (printer trace)
  (let ((success nil)
        (total-time 0)
        prev-numberify)
    (unwind-protect
         (loop with count = (length *solvers*)
               for i upfrom 1
               for (solver numberify) in *solvers* do
           (multiple-value-bind (obj status values time)
               (%solve-with-solver solver
                                   (and (not (eq numberify prev-numberify))
                                        (lambda (stream)
                                          (funcall printer numberify stream)))
                                   success)
             (setf prev-numberify numberify)
             (format trace "~A: ~,3F sec ~F (~(~A~))~%"
                     solver time (double obj) status)
             ;; FIXME: are basis file sane when !optimal/unbounded?
             (setf success (member status '(:optimal :infeasible
                                            :unbounded)))
             (incf total-time time)
             (when (= i count)
               (return (values status obj values total-time)))))
      (flet ((try-delete (x)
               (when (probe-file x)
                 (delete-file x))))
        (try-delete (format nil "~A.lp"
                            *instance-stem*))
        (try-delete (format nil "~A.bas"
                            *instance-stem*))))))

(defun solve-with-printer (printer &key (trace *trace-output*)
                                     inexact double-only)
  (let ((*solvers* (cond (double-only
                          (subseq *solvers* 0 1))
                         (inexact
                          (butlast *solvers*))
                         (t *solvers*))))
    (%solve printer trace)))

(defun solve (&key model (trace *trace-output*)
                inexact double-only)
  (setf model (or model *model*))
  (multiple-value-bind (status obj values time)
      (solve-with-printer (lambda (numberify stream)
                            (%print-model model numberify stream))
                          :trace trace
                          :inexact inexact
                          :double-only double-only)
    (declare (ignore obj))
    (let ((table  (make-hash-table))
          (obj    nil)
          (id-var (id-var model))
          (*read-default-float-format* 'double-float))
      (maphash (lambda (name value)
                 (cond ((equal name "obj_var")
                        (assert (not obj))
                        (setf obj value))
                       (t
                        (let* ((idx (parse-integer name :start 1))
                               (var (aref id-var idx)))
                          (assert (not (gethash var table)))
                          (setf (gethash var table) value)))))
               values)
      (values status obj table time))))

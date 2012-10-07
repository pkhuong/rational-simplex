(defvar *float-mode* 'double-float)

(defun %double-float-bits (x)
  (let* ((x    (float x 1d0))
         (bits (logior (ash (sb-kernel:double-float-high-bits x) 32)
                       (sb-kernel:double-float-low-bits x))))
    ;; convert from sign-magnitude
    (if (logbitp 63 bits)
        (- (ldb (byte 63 0) bits))
        bits)))

(defun %double-float-from-bits (x &optional delta)
  ;; revert into sign/magnitude
  (when delta (incf x delta))
  (when (minusp x)
    (setf x (logior (ash -1 63)
                    (ldb (byte 63 0) (- x)))))
  (sb-kernel:make-double-float (ash x -32)
                               (ldb (byte 32 0) x)))


(defun float-bits (x)
  (ecase *float-mode*
    (double-float
     (%double-float-bits x))
    (single-float
     (let ((bits (sb-kernel:single-float-bits (float x 1s0))))
       (if (logbitp 31 bits)
           (- (ldb (byte 31 0) bits))
           bits)))))

(defun float-from-bits (x &optional delta)
  (ecase *float-mode*
    (double-float
     (%double-float-from-bits x delta))
    (single-float
     (when delta
       (incf x delta))
     (when (minusp x)
       (setf x (logior (ash -1 31)
                       (ldb (byte 63 0) (- x)))))
     (sb-kernel:make-single-float x))))

(defun round-to-float (x)
  (rational (coerce x *float-mode*)))

(defun floatify (x)
  (coerce x *float-mode*))

(defun split-float (x)
  (let ((dx (round-to-float x)))
    (cond ((= x dx)
           nil)
          ((< x dx)
           (values (float-from-bits (1- (float-bits dx)))
                   dx))
          (t
           (values dx
                   (float-from-bits (1+ (float-bits dx))))))))

(defvar *precision* 64)
(defvar *max-accuracy* 4096)

(defun pull-bits (x &optional (*precision* *precision*)
                  &aux (x (if (floatp x)
                              (rational x)
                              x)))
  (labels ((done (a shift)
             (let* ((len (integer-length a))
                    (delta (- len *precision*)))
               (assert (>= len *precision*))
               (return-from pull-bits
                 (/ (round a (ash 1 delta))
                    (ash 1 (- shift delta))))))
           (pull-bits (shift)
             (let ((a (computable-reals:approx-r x shift)))
               (cond ((zerop a))
                     ((> (integer-length a) *precision*)
                      (done a shift))
                     (t
                      (pull-bits (+ shift *precision*)))))))
    (dotimes (i (integer-length *max-accuracy*) 0)
      (pull-bits (ash 1 i)))))

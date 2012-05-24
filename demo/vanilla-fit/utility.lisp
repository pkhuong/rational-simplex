(defun double-float-bits (x)
  (let* ((x    (float x 1d0))
         (bits (logior (ash (sb-kernel:double-float-high-bits x) 32)
                       (sb-kernel:double-float-low-bits x))))
    ;; convert from sign-magnitude
    (if (logbitp 63 bits)
        (- (ldb (byte 63 0) bits))
        bits)))

(defun double-float-from-bits (x)
  ;; revert into sign/magnitude
  (when (minusp x)
    (setf x (logior (ash -1 63)
                    (ldb (byte 63 0) (- x)))))
  (sb-kernel:make-double-float (ash x -32)
                               (ldb (byte 32 0) x)))

(defun round-to-double (x)
  (rational (float x 1d0)))

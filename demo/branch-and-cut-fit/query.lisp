CL-USER> (defparameter *worker*
           (sb-thread:make-thread
            (lambda ()
              (let ((*standard-output* sb-impl::*stdout*)
                    (*trace-output*    sb-impl::*stdout*))
                #+nil
                (print-hash-table
                 (time 
                  (enumerate-pareto-front 
                   16 -2 2
                   #'computable-reals:sin-r
                   #'computable-reals:cos-r
                   (lambda (x)
                     (computable-reals:-r (computable-reals:sin-r x)))))
                 "/tmp/sin-16")
                (print-hash-table
                 (time 
                  (enumerate-pareto-front 
                   16 0 1
                   (lambda (x)
                     (computable-reals:log-r (computable-reals:+r x 1)))
                   (lambda (x)
                     (let ((x (computable-reals:+r x 1)))
                       (computable-reals:/R x)))
                   (lambda (x)
                     (let ((x (computable-reals:+r x 1)))
                       (computable-reals:-r 
                        (computable-reals:/R
                         (computable-reals:*r x x)))))))
                 "/tmp/log1px-16")
                #+nil
                (print-hash-table
                 (time 
                  (enumerate-pareto-front 
                   16 1 2
                   (lambda (x)
                     (computable-reals:log-r x))
                   (lambda (x)
                     (computable-reals:/R x))
                   (lambda (x)
                     (computable-reals:-r 
                      (computable-reals:/R
                       (computable-reals:*r x x))))))
                 "/tmp/log-16")
                #+nil
                (print-hash-table
                 (time 
                  (enumerate-pareto-front 
                   16 0 1
                   #'computable-reals:exp-r
                   #'computable-reals:exp-r
                   #'computable-reals:exp-r))
                 "/tmp/exp-16")))))

*WORKER*

CL-USER> (defparameter *worker2*
           (sb-thread:make-thread
            (lambda ()
              (sb-thread:join-thread *worker*)
              (let ((*standard-output* sb-impl::*stdout*)
                    (*trace-output*    sb-impl::*stdout*))
                (print-hash-table
                 (time 
                  (enumerate-pareto-front 
                   16
                   (- (double-float-from-bits (double-float-bits (/ pi 2)) 2))
                   (double-float-from-bits (double-float-bits (/ pi 2)) 2)
                   #'computable-reals:sin-r
                   #'computable-reals:cos-r
                   (lambda (x)
                     (computable-reals:-r (computable-reals:sin-r x)))))
                 "/tmp/sin-narrow-16")
                (print-hash-table
                 (time 
                  (enumerate-pareto-front 
                   16 -1 1
                   #'computable-reals:exp-r
                   #'computable-reals:exp-r
                   #'computable-reals:exp-r))
                 "/tmp/exp-wide-16")))))
*WORKER2*
CL-USER> (defparameter *worker*
           (sb-thread:make-thread
            (lambda ()
              (let ((*standard-output* sb-impl::*stdout*)
                    (*trace-output*    sb-impl::*stdout*))
                (print-hash-table
                 (time 
                  (enumerate-pareto-front 
                   16 -1 1
                   #'computable-reals:atan-r
                   (lambda (x)
                     (computable-reals:/r (computable-reals:+r 
                                           1
                                           (computable-reals:*r x x))))
                   (lambda (x)
                     (let ((x^2+1 (computable-reals:+r
                                   1 (computable-reals:*r x x))))
                       (computable-reals:-r
                        (computable-reals:/r (computable-reals:*r 2 x)
                                             x^2+1))))))                   
                 "/tmp/arctan-16")))))
*WORKER*
CL-USER> (defparameter *worker*
           (sb-thread:make-thread
            (lambda ()
              (let ((*standard-output* sb-impl::*stdout*)
                    (*trace-output*    sb-impl::*stdout*))
                (print-hash-table
                 (time 
                  (enumerate-pareto-front 
                   16 0 1
                   #'computable-reals:atan-r
                   (lambda (x)
                     (computable-reals:/r (computable-reals:+r 
                                           1
                                           (computable-reals:*r x x))))
                   (lambda (x)
                     (let ((x^2+1 (computable-reals:+r
                                   1 (computable-reals:*r x x))))
                       (computable-reals:-r
                        (computable-reals:/r (computable-reals:*r 2 x)
                                             x^2+1))))))                   
                 "/tmp/arctan-2-16")))))
*WORKER*
CL-USER> (defparameter *worker*
           (sb-thread:make-thread
            (lambda ()
              (let ((*standard-output* sb-impl::*stdout*)
                    (*trace-output*    sb-impl::*stdout*))
                (print-hash-table
                 (time 
                  (enumerate-pareto-front 
                   16 0 1
                   (lambda (x)
                     (let ((x (computable-reals:+r 1 x)))
                       (computable-reals:log-r x 2)))
                   (lambda (x)
                     (let ((x (computable-reals:+r 1 x)))
                       (computable-reals:/r
                        (computable-reals:*r (computable-reals:log-r 2)
                                             x))))
                   (lambda (x)
                     (let ((x (computable-reals:+r 1 x)))
                       (computable-reals:-r
                        (computable-reals:/r
                         (computable-reals:*r (computable-reals:log-r 2)
                                              x x)))))))
                 "/tmp/lg1px-16")))))
*WORKER*
CL-USER> (Setf *float-mode* 'single-float)
SINGLE-FLOAT
CL-USER> (defparameter *worker*
           (sb-thread:make-thread
            (lambda ()
              (let ((*standard-output* sb-impl::*stdout*)
                    (*trace-output*    sb-impl::*stdout*))
                (print-hash-table
                 (time 
                  (enumerate-pareto-front 
                   16 0 1
                   (lambda (x)
                     (computable-reals:log-r (computable-reals:+r x 1)))
                   (lambda (x)
                     (let ((x (computable-reals:+r x 1)))
                       (computable-reals:/R x)))
                   (lambda (x)
                     (let ((x (computable-reals:+r x 1)))
                       (computable-reals:-r 
                        (computable-reals:/R
                         (computable-reals:*r x x)))))))
                 "/tmp/single/log1px-16")
                (print-hash-table
                 (time 
                  (enumerate-pareto-front 
                   16 1 2
                   (lambda (x)
                     (computable-reals:log-r x))
                   (lambda (x)
                     (computable-reals:/R x))
                   (lambda (x)
                     (computable-reals:-r 
                      (computable-reals:/R
                       (computable-reals:*r x x))))))
                 "/tmp/single/log1px-16")
                (print-hash-table
                 (time 
                  (enumerate-pareto-front 
                   16
                   (- (float-from-bits (float-bits (/ pi 2)) 2))
                   (float-from-bits (float-bits (/ pi 2)) 2)
                   #'computable-reals:sin-r
                   #'computable-reals:cos-r
                   (lambda (x)
                     (computable-reals:-r (computable-reals:sin-r x)))))
                 "/tmp/single/sin-narrow-16")
                (print-hash-table
                 (time 
                  (enumerate-pareto-front 
                   16
                   (- (float-from-bits (float-bits (/ pi 2)) 2))
                   (float-from-bits (float-bits (/ pi 2)) 2)
                   #'computable-reals:cos-r
                   (lambda (x)
                     (computable-reals:-r (computable-reals:sin-r x)))
                   (lambda (x)
                     (computable-reals:-r (computable-reals:cos-r x)))))
                 "/tmp/single/cos-narrow-16")
                (print-hash-table
                 (time 
                  (enumerate-pareto-front 
                   16 -1 1
                   #'computable-reals:exp-r
                   #'computable-reals:exp-r
                   #'computable-reals:exp-r))
                 "/tmp/single/exp-wide-16")
                (print-hash-table
                 (time 
                  (enumerate-pareto-front 
                   16 0 1
                   #'computable-reals:exp-r
                   #'computable-reals:exp-r
                   #'computable-reals:exp-r))
                 "/tmp/single/exp-16")
                (print-hash-table
                 (time 
                  (enumerate-pareto-front 
                   16 0 1
                   #'computable-reals:atan-r
                   (lambda (x)
                     (computable-reals:/r (computable-reals:+r 
                                           1
                                           (computable-reals:*r x x))))
                   (lambda (x)
                     (let ((x^2+1 (computable-reals:+r
                                   1 (computable-reals:*r x x))))
                       (computable-reals:-r
                        (computable-reals:/r (computable-reals:*r 2 x)
                                             x^2+1))))))                   
                 "/tmp/single/arctan-16")
                (print-hash-table
                 (time 
                  (enumerate-pareto-front 
                   16 -1 1
                   #'computable-reals:atan-r
                   (lambda (x)
                     (computable-reals:/r (computable-reals:+r 
                                           1
                                           (computable-reals:*r x x))))
                   (lambda (x)
                     (let ((x^2+1 (computable-reals:+r
                                   1 (computable-reals:*r x x))))
                       (computable-reals:-r
                        (computable-reals:/r (computable-reals:*r 2 x)
                                             x^2+1))))))                   
                 "/tmp/single/arctan-wide-16")

                (print-hash-table
                 (time 
                  (enumerate-pareto-front 
                   16 0 1
                   (lambda (x)
                     (let ((x (computable-reals:+r 1 x)))
                       (computable-reals:log-r x 2)))
                   (lambda (x)
                     (let ((x (computable-reals:+r 1 x)))
                       (computable-reals:/r
                        (computable-reals:*r (computable-reals:log-r 2)
                                             x))))
                   (lambda (x)
                     (let ((x (computable-reals:+r 1 x)))
                       (computable-reals:-r
                        (computable-reals:/r
                         (computable-reals:*r (computable-reals:log-r 2)
                                              x x)))))))
                 "/tmp/single/lg1px-16")))))
*WORKER*

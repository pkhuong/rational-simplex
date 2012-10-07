(asdf:defsystem "bac-minimax-polynomial"
  :version "0.0.1"
  :licence "BSD"
  :description "Branch and cut solver for minimax polynomial approximation"
  :serial t
  :depends-on ("rational-simplex" "computable-reals" "sb-concurrency"
                                  "sb-md5")
  :components ((:file "utility")
               (:file "linf-fit")
               (:file "newton")
               (:file "find-extrema")
               (:file "driver")
               (:file "print-pareto")))

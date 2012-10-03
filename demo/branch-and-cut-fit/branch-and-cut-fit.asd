(asdf:defsystem "bac-minimax-polynomial"
  :version "0.0.1"
  :licence "BSD"
  :description "Branch and cut solver for minimax polynomial approximation"
  :serial t
  :depends-on ("rational-simplex" "computable-reals")
  :components ((:file "utility")
               (:file "linf-fit")
               (:file "newton")
               (:file "find-extrema")
               (:file "driver")))

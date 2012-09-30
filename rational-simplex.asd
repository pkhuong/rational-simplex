(asdf:defsystem "rational-simplex"
  :version "0.0.1"
  :licence "BSD"
  :description "Tiny modeling language and wrapper around QSopt-Exact"
  :serial t
  :depends-on ("cl-ppcre" #-sbcl "trivial-shell" "alexandria")
  :components ((:file "rational-simplex")))

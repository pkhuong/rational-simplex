RATIONAL-SIMPLEX
================

This one-file system provides two things: a tiny modeling framework
for linear programs, and a wrapper around Daniel Espinoza et al's
QSopt-Exact, an exact (with rational arithmetic) solver for linear
optimization problems.

Installation
------------

The solver requires an installation of
[QSopt-Exact](http://www.dii.uchile.cl/~daespino/ESolver_doc/main.html).
The program depends on
[EGlib](http://www.dii.uchile.cl/~daespino/EGlib_doc/main.html), an
utility library, which has a bug that results in erroneous variable
value output for very large fractions or floats.  Both programs are
GPL.

If you intend to work with large denominators, you'll need to build
both programs from sources (it only takes a few minutes).  Applying
`EGlib-long-lines.patch` with `patch -p0 < EGlib-long-lines.patch`
from `EGlib-2.6.20/` before building EGlib (and QSopt_ex) should
work... It's an ugly workaround, you probably don't want to look too
closely (:

`*solvers-path*` should point to the directory holding the QSopt_ex
executables.  It defaults to
`[directory from which simplex.lisp is loaded]/bin`, which is probably
not what you want if you're loading it with asdf or quicklisp.

`*instance-stem*` defaults to `"/tmp/rational-simplex-instance"`;
temporary files `/tmp/rational-simplex-instance.lp` and
`/tmp/rational-simplex-instance.bas` will be created (and
overwritten!).  You probably want to change that to something more
unique (and private) if you work on a shared machine.

Other dependencies: `cl-ppcre` and `trivial-shell` on `#-sbcl`
platforms (in which case, beware spaces in paths...).

Modeling language
-----------------

The base object is the `model`.  It stores the current objective,
sense (should the objective function be minimised or maximised),
variables and constraints.  `with-model (&key name sense)` should be
used in most cases (`name` is a string or nil, and `sense` is
`:minimize`, the default, or `:maximise`).

In the scope of `with-model`, `var &key name lower upper obj`
instantiates a new variables associated with the current model.  The
lower bound defaults to 0, and the upper bound to none (`nil`).  `obj`
is the variable's coefficient in the objective function, and defaults
to 0.

Variables and real values can be composed together into immutable
`linear-expression`s.  `linexpr &optional constant &rest {coef var}*`
creates a fresh linear expression corresponding to `constant +
coef1*var1 + ...`.  `add` or `sub` can be used to add or subtract
reals, variables, or linear expressions; the optional third argument
specifies a scaling value for the second argument.  `scale` will scale
it argument by a real, while `remove-var` can be used to completely
remove a variable from a linear expression. `addf`, `subf`, `scalef`,
`remove-varf` are convenience modify macros around the functions.

The function `constraint lhs cmp rhs` returns a new constraints that
asserts a relationship between two linear expressions (or variables or
reals).  `cmp` can be `<=`, `>=` or `=`.  `add-constraint` and
`del-constraint` add/delete a constraint to/from a model.  `constrain
lhs cmp rhs &optional model` or its synonym `post` create a constraint
and add it to the model (defaults to `*model*`) in a single call.  The
return value is the new constraint itself, to easily `del-constraint`
it.

`constraints` is a place that ensures that a model always has unique
ownership over the constraint store (copies are stored/returned).
This can be used to simplify branch-and-bound type algorithms.

`print-float-model` and `print-rational-model` can be used to print
the current model in lp format.  Values in the model, linear
expressions, etc. are converted to rationals as early as possible, so
`print-rational-model` does not lose precision.  Most solvers work
with floating point values; in that case `print-float-model` is more
appropriate.

Solver
------

`solve &key model trace inexact double-only` is used to solve a model
specified in the framework described above.  `model` defaults to
`*model*`.  `trace` specifies where short (a few lines line per solve)
logging should be printed and defaults to `*trace-output*`.  `inexact`
defaults to `nil`; if true, the final solve in rational is skipped.
`double-only` defaults to `nil`; if true, only the initial solve with
double arithmetic is performed.  The return values are a status
(`:optimal`, `:infeasible`, `:unbounded` or `unsolved`), the objective
value, a hash table from variables to values, and the total solution
time.

This is a wrapper around `solve-with-printer printer &key trace
inexact double-only`.  Instead of a model, this function receives a
printer function which, givena function with which to canonicalize
numbers (e.g. one that converts any real to a float) and a stream,
prints an lp model to the stream.  The return values are the last
executed solver's status, the reported objective value, a hash table
mapping variable names to values, and the total solution time.  Note
that the objective value is reported with little precision;
associating the expression with a bogus variable and optimising on
that variable will give full precision in the objective value.

Example usage
-------------

This solves a silly 3-variable linear program, and reports the
objective value and decision variables' values at an optimal
solution.

    CL-USER> (lp:with-model (:name "small" :sense :maximize)
               (let ((x (lp:var :name "x" :lower 2 :obj 3))
                     (y (lp:var :name "y" :lower nil :obj 2))
                     (z (lp:var :name "z" :lower 1 :upper 10 :obj 4)))
                 (lp:post (lp:linexpr 0 3 x 2 y 1 z) '<= 12)
                 (lp:post (lp:add y x 5) '<= 10)
                 (multiple-value-bind (status obj values)
                     (lp:solve :double-only nil :inexact nil)
                   (assert (eql status :optimal))
                   (values obj (mapcar (lambda (var)
                                         (cons var (gethash var values)))
                                       (list x y z))))))
    dbl_solver: 0.000 sec 0.0 (optimal)
    ldbl_solver: 0.000 sec 0.0 (optimal)
    float128_solver: 0.000 sec 0.0 (optimal)
    mpf_solver: 0.010 sec 0.0 (optimal)
    mpq_solver: 0.010 sec 0.0 (optimal)
    42
    ((#<RATIONAL-SIMPLEX:VAR x [2.0:+inf]> . 2)
     (#<RATIONAL-SIMPLEX:VAR y [-inf:+inf]> . -2)
     (#<RATIONAL-SIMPLEX:VAR z [1.0:10.0]> . 10))

Implementation
--------------

QSopt-Exact comes with a full adaptive-precision MIP solver... but it
solution output routines seem broken on my linux, and I didn't
particularly feel like going down that rabbit hole.  Instead, I call
LP solvers of increasing precision in sequence: double, long doubles,
128 bit quads, multi-precision floats and finally rationals.  The
reason this works (and is interesting) is that the simplex is a
combinatorial algorithm: the correctness of each step only depends on
the *comparison* between two values (i.e. is `x >= y`).  Better: the
state determined at each step isn't numeric, but rather a basis (a set
of variable or constraints).  Thus, we can save the basis found by the
double solver, and use it to warm start the more precise long double
solver, etc.  Our hope is that the float-based solvers will converge
to the optimal solution, or very nearly, leaving only a few iterations
until the final exact, slow, solver itself declares victory.

;;; -*- Mode: lisp -*-

(defpackage #:quad-double-internal
  (:use #:cl #+cmu #:extensions)
  (:nicknames #:qdi)
  (:export #:%quad-double
	   #:read-qd
	   #:add-qd
	   #:add-qd-d
	   #:cmu #:add-qd-dd
	   #:add-d-qd
	   #:sub-qd
	   #:sub-qd-d
	   #:cmu #:sub-qd-dd
	   #:sub-d-qd
	   #:neg-qd
	   #:mul-qd
	   #:mul-qd-d
	   #:sqr-qd
	   #:div-qd
	   #:div-qd-d
	   #+cmu #:div-qd-dd
	   #:make-qd-d
	   #+cmu #:make-qd-dd
	   #:integer-decode-qd
	   #:npow
	   #:qd-0
	   #:qd-1
	   #:qd-2
	   #:qd-3
	   #:qd-parts
	   #:+qd-one+
	   #:+qd-zero+
	   #:+qd-pi+
	   ;; Functions
	   #:abs-qd
	   #:sqrt-qd
	   #:log-qd
	   #:log1p-qd
	   #:exp-qd
	   #:sin-qd
	   #:cos-qd
	   #:tan-qd
	   #:sincos-qd
	   #:asin-qd
	   #:acos-qd
	   #:atan-qd
	   #:atan2-qd
	   #:sinh-qd
	   #:cosh-qd
	   #:tanh-qd
	   #:asinh-qd
	   #:acosh-qd
	   #:atanh-qd
	   #:qd-=
	   #:qd->
	   #:qd-<
	   #:qd->=
	   #:qd-<=
	   #:zerop-qd
	   #:plusp-qd
	   #:minusp-qd
	   #:integer-decode-qd
	   #:decode-float-qd
	   #:scale-float-qd
	   #:ffloor-qd
	   )
  #+cmu
  (:import-from #:c
		#:two-sum
		#:quick-two-sum
		#:two-prod
		#:two-sqr))

(defpackage #:quad-double
  (:use #:cl #:quad-double-internal)
  (:nicknames #:qd)
  (:shadow #:+
	   #:-
	   #:*
	   #:/
	   #:1+
	   #:1-
	   #:zerop
	   #:plusp
	   #:minusp
	   #:abs
	   #:sqrt
	   #:log
	   #:exp
	   #:sin
	   #:cos
	   #:tan
	   #:asin
	   #:acos
	   #:atan
	   #:sinh
	   #:cosh
	   #:tanh
	   #:asinh
	   #:acosh
	   #:atanh
	   #:expt
	   #:=
	   #:/=
	   #:<
	   #:>
	   #:<=
	   #:>=
	   #:complex
	   #:integer-decode-float
	   #:decode-float
	   #:scale-float
	   #:float
	   #:floor
	   #:ffloor
	   #:ceiling
	   #:fceiling
	   #:truncate
	   #:ftruncate
	   #:round
	   #:fround
	   #:realpart
	   #:imagpart
	   #:conjugate
	   #:float-sign
	   #:qd-format-exp
	   #:max
	   #:min
	   #:cis
	   #:phase
	   #:signum
	   )
  (:export #:+
	   #:-
	   #:*
	   #:/
	   #:1+
	   #:1-
	   #:zerop
	   #:plusp
	   #:minusp
	   #:abs
	   #:sqrt
	   #:log
	   #:exp
	   #:sin
	   #:cos
	   #:tan
	   #:asin
	   #:acos
	   #:atan
	   #:sinh
	   #:cosh
	   #:tanh
	   #:asinh
	   #:acosh
	   #:atanh
	   #:expt
	   #:=
	   #:/=
	   #:<
	   #:>
	   #:<=
	   #:>=
	   #:complex
	   #:integer-decode-float
	   #:decode-float
	   #:scale-float
	   #:float
	   #:floor
	   #:ffloor
	   #:ceiling
	   #:fceiling
	   #:truncate
	   #:ftruncate
	   #:round
	   #:fround
	   #:realpart
	   #:imagpart
	   #:conjugate
	   #:float-sign
	   #:qd-format-exp
	   #:max
	   #:min
	   #:cis
	   #:phase
	   #:signum
	   )
  ;; Constants
  (:export #:+pi+)
  ;; CMUCL supports infinities.
  #+cmu
  (:export #:+quad-double-float-positive-infinity+
	   #:+quad-double-float-negative-infinity+))

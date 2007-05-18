;;; -*- Mode: lisp -*-

(defpackage "QUAD-DOUBLE-INTERNAL"
  (:use "CL" "EXTENSIONS")
  (:nicknames "QDI")
  (:export "%QUAD-DOUBLE"
	   "READ-QD"
	   "ADD-QD"
	   "ADD-QD-D"
	   "ADD-D-QD"
	   "SUB-QD"
	   "SUB-QD-D"
	   "SUB-D-QD"
	   "NEG-QD"
	   "MUL-QD"
	   "MUL-QD-D"
	   "DIV-QD"
	   "DIV-QD-D"
	   "MAKE-QD-D"
	   "MAKE-QD-DD"
	   "INTEGER-DECODE-QD"
	   "NPOW"
	   "QD-0"
	   "+QD-ONE+"
	   "+QD-ZERO+"
	   "+QD-PI+"
	   ;; Functions
	   "SQRT-QD"
	   "LOG-QD"
	   "EXP-QD"
	   "SIN-QD"
	   "COS-QD"
	   "TAN-QD"
	   "ASIN-QD"
	   "ACOS-QD"
	   "ATAN-QD"
	   "ATAN2-QD"
	   "SINH-QD"
	   "COSH-QD"
	   "TANH-QD"
	   "ASINH-QD"
	   "ACOSH-QD"
	   "ATANH-QD"
	   "QD-="
	   "QD->"
	   "QD-<"
	   "QD->="
	   "QD-<="
	   "ZEROP-QD"
	   "PLUSP-QD"
	   "MINUSP-QD"
	   "INTEGER-DECODE-QD"
	   "DECODE-FLOAT-QD"
	   "SCALE-FLOAT-QD"
	   ))

(defpackage "QUAD-DOUBLE"
  (:use "CL" "QUAD-DOUBLE-INTERNAL")
  (:nicknames "QD")
  (:shadow "+"
	   "-"
	   "*"
	   "/"
	   "1+"
	   "1-"
	   "ZEROP"
	   "PLUSP"
	   "MINUSP"
	   "SQRT"
	   "LOG"
	   "EXP"
	   "SIN"
	   "COS"
	   "TAN"
	   "ASIN"
	   "ACOS"
	   "ATAN"
	   "SINH"
	   "COSH"
	   "TANH"
	   "ASINH"
	   "ACOSH"
	   "ATANH"
	   "EXPT"
	   "="
	   "/="
	   "<"
	   ">"
	   "<="
	   ">="
	   "COMPLEX"
	   "INTEGER-DECODE-FLOAT"
	   "DECODE-FLOAT"
	   "SCALE-FLOAT"
	   ))

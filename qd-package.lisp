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
	   "+QD-ZERO+"))

(defpackage "QUAD-DOUBLE"
  (:use "CL" "QUAD-DOUBLE-INTERNAL")
  (:nicknames "QD")
  (:shadow "+"
	   "-"
	   "*"
	   "/"))

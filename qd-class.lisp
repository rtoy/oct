(in-package "QD")

(define-symbol-macro * cl:*)
(define-symbol-macro - cl:-)
(define-symbol-macro / cl:/)

(defclass quad-double ()
  ((qd :initform +qd-zero+
       :reader qd-value
       :initarg :value
       :type %quad-double)))

(defmethod print-object ((qd quad-double) stream)
  (format stream "#q~/qdi::qd-format/" (qd-value qd)))

(defmethod make-qd ((x real))
  (make-instance 'quad-double :value (make-qd-d (float x 1d0))))

(defmethod make-qd ((x quad-double))
  (make-instance 'quad-double :value (qd-value x)))
  
(defun qd-class-reader (stream subchar arg)
  (declare (ignore subchar arg))
  (make-instance 'quad-double :value (read-qd stream)))

;; Yow!  We redefine the #q reader that is in qd-io.lisp to read in
;; and make a real quad-double float, instead of the hackish
;; %quad-double.
(set-dispatch-macro-character #\# #\Q #'qd-class-reader)

(defmethod add1 ((a number))
  (cl::1+ a))

(defmethod add1 ((a quad-double))
  (make-instance 'quad-double :value (add-qd-d (qd-value a) 1d0)))

(defmethod sub1 ((a number))
  (cl::1- a))

(defmethod sub1 ((a quad-double))
  (make-instance 'quad-double :value (sub-qd-d (qd-value a) 1d0)))

(declaim (inline 1+ 1-))

(defun 1+ (x)
  (add1 x))

(defun 1- (x)
  (sub1 x))

(defmethod two-arg-+ ((a quad-double) (b quad-double))
  (make-instance 'quad-double :value (add-qd (qd-value a) (qd-value b))))

(defmethod two-arg-+ ((a quad-double) (b real))
  (make-instance 'quad-double :value (add-qd-d (qd-value a) (float b 1d0))))

(defmethod two-arg-+ ((a real) (b quad-double))
  (make-instance 'quad-double :value (add-d-qd (float a 1d0) (qd-value b))))

(defmethod two-arg-+ ((a number) (b number))
  (cl:+ a b))

(defun + (&rest args)
  (if (null args)
      0
      (do ((args (cdr args) (cdr args))
	   (res (car args)
		(two-arg-+ res (car args))))
	  ((null args) res))))

(defmethod two-arg-- ((a quad-double) (b quad-double))
  (make-instance 'quad-double :value (sub-qd (qd-value a) (qd-value b))))

(defmethod two-arg-- ((a quad-double) (b real))
  (make-instance 'quad-double :value (sub-qd-d (qd-value a) (float b 1d0))))

(defmethod two-arg-- ((a real) (b quad-double))
  (make-instance 'quad-double :value (sub-d-qd (float a 1d0) (qd-value b))))

(defmethod two-arg-- ((a number) (b number))
  (cl:- a b))

(defmethod unary-minus ((a number))
  (cl:- a))

(defmethod unary-minus ((a quad-double))
  (make-instance 'quad-double :value (neg-qd (qd-value a))))

(defun - (number &rest more-numbers)
  (if more-numbers
      (do ((nlist more-numbers (cdr nlist))
	   (result number))
	  ((atom nlist) result)
         (declare (list nlist))
	 (setq result (two-arg-- result (car nlist))))
      (unary-minus number)))


(defmethod two-arg-* ((a quad-double) (b quad-double))
  (make-instance 'quad-double :value (mul-qd (qd-value a) (qd-value b))))

(defmethod two-arg-* ((a quad-double) (b real))
  (make-instance 'quad-double :value (mul-qd-d (qd-value a) (float b 1d0))))

(defmethod two-arg-* ((a real) (b quad-double))
  (make-instance 'quad-double :value (mul-qd-d (qd-value b) (float a 1d0))))

(defmethod two-arg-* ((a number) (b number))
  (cl:* a b))

(defun * (&rest args)
  (if (null args)
      1
      (do ((args (cdr args) (cdr args))
	   (res (car args)
		(two-arg-* res (car args))))
	  ((null args) res))))

(defmethod two-arg-/ ((a quad-double) (b quad-double))
  (make-instance 'quad-double :value (div-qd (qd-value a) (qd-value b))))

(defmethod two-arg-/ ((a quad-double) (b real))
  (make-instance 'quad-double :value (div-qd-d (qd-value a) (float b 1d0))))

(defmethod two-arg-/ ((a real) (b quad-double))
  (make-instance 'quad-double :value (div-qd (make-qd-d (float a 1d0))
					     (qd-value b))))

(defmethod two-arg-/ ((a number) (b number))
  (cl:/ a b))

(defmethod unary-divide ((a number))
  (cl:/ a))

(defmethod unary-divide ((a quad-double))
  (make-instance 'quad-double :value (div-qd +qd-one+ (qd-value a))))

(defun / (number &rest more-numbers)
  (if more-numbers
      (do ((nlist more-numbers (cdr nlist))
	   (result number))
	  ((atom nlist) result)
         (declare (list nlist))
	 (setq result (two-arg-/ result (car nlist))))
      (unary-divide number)))

(defmethod qlog ((a real) &optional b)
  (if b
      (cl:log a b)
      (cl:log a)))

(defmethod qlog ((a quad-double) &optional b)
  (make-instance 'quad-double
		 :value (if b
			    (/ (log-qd (qd-value a))
			       (if (realp b)
				   (cl:log b)
				   (log-qd (make-qd-d (float b 1d0)))))
			    (log-qd (qd-value a)))))

(declaim (inline log))
(defun log (a &optional b)
  (qlog a b))

(macrolet ((frob (name)
	     (let ((method-name (intern (concatenate 'string "Q" (symbol-name name))))
		   (cl-name (intern (symbol-name name) :cl))
		   (qd-name (intern (concatenate 'string (symbol-name name) "-QD"))))
	       `(progn
		 (defmethod ,method-name ((x number))
		   (,cl-name x))
		 (defmethod ,method-name ((x quad-double))
		   (make-instance 'quad-double :value (,qd-name (qd-value x))))
		 (declaim (inline ,name))
		 (defun ,name (x)
		   (,method-name x))))))
  (frob sqrt)
  (frob exp)
  (frob sin)
  (frob cos)
  (frob tan)
  (frob asin)
  (frob acos)
  (frob sinh)
  (frob cosh)
  (frob tanh)
  (frob asinh)
  (frob acosh)
  (frob atanh))

(defmethod qatan ((y real) &optional x)
  (cl:atan y x))

(defmethod qatan ((y quad-double) &optional x)
  (make-instance 'quad-double
		 :value
		 (if x
		     (atan2-qd (qd-value y) (qd-value x))
		     (atan-qd (qd-value y)))))

(defun atan (y &optional x)
  (qatan y x))



(macrolet
    ((frob (op)
       (let ((method (intern (concatenate 'string "TWO-ARG-" (symbol-name op))))
	     (cl-fun (find-symbol (symbol-name op) :cl))
	     (qd-fun (intern (concatenate 'string "QD-" (symbol-name op))
			     (find-package :qdi))))
	 `(progn
	    (defmethod ,method ((a real) (b real))
	      (,cl-fun a b))
	    (defmethod ,method ((a quad-double) (b real))
	      (,qd-fun (qd-value a) (make-qd-d (float b 1d0))))
	    (defmethod ,method ((a real) (b quad-double))
	      (,qd-fun (make-qd-d (float a 1d0)) (qd-value b)))))))
  (frob <)
  (frob >)
  (frob <=)
  (frob >=))

(defmethod two-arg-= ((a number) (b number))
  (cl:= a b))
(defmethod two-arg-= ((a quad-double) (b number))
  (if (realp b)
      (qdi::qd-= (qd-value a) (make-qd-d (float b 1d0)))
      nil))
(defmethod two-arg-= ((a number) (b quad-double))
  (if (realp a)
      (qdi::qd-= (make-qd-d (float a 1d0)) (qd-value b))
      nil))

(defun = (number &rest more-numbers)
  "Returns T if all of its arguments are numerically equal, NIL otherwise."
  (declare (optimize (safety 2))
	   (dynamic-extent more-numbers))
  (do ((nlist more-numbers (cdr nlist)))
      ((atom nlist) T)
    (declare (list nlist))
    (if (not (two-arg-= (car nlist) number))
	(return nil))))

(defun /= (number &rest more-numbers)
  "Returns T if no two of its arguments are numerically equal, NIL otherwise."
  (declare (optimize (safety 2))
	   (dynamic-extent more-numbers))
  (do* ((head number (car nlist))
	(nlist more-numbers (cdr nlist)))
       ((atom nlist) t)
    (declare (list nlist))
    (unless (do* ((nl nlist (cdr nl)))
		 ((atom nl) T)
	      (declare (list nl))
	      (if (two-arg-= head (car nl))
		  (return nil)))
      (return nil))))

(defun < (number &rest more-numbers)
  "Returns T if its arguments are in strictly increasing order, NIL otherwise."
  (declare (optimize (safety 2))
	   (dynamic-extent more-numbers))
  (do* ((n number (car nlist))
	(nlist more-numbers (cdr nlist)))
       ((atom nlist) t)
     (declare (list nlist))
     (if (not (two-arg-< n (car nlist))) (return nil))))

(defun > (number &rest more-numbers)
  "Returns T if its arguments are in strictly decreasing order, NIL otherwise."
  (declare (optimize (safety 2))
	   (dynamic-extent more-numbers))
  (do* ((n number (car nlist))
	(nlist more-numbers (cdr nlist)))
       ((atom nlist) t)
     (declare (list nlist))
     (if (not (two-arg-> n (car nlist))) (return nil))))

(defun <= (number &rest more-numbers)
  "Returns T if arguments are in strictly non-decreasing order, NIL otherwise."
  (declare (optimize (safety 2))
	   (dynamic-extent more-numbers))
  (do* ((n number (car nlist))
	(nlist more-numbers (cdr nlist)))
       ((atom nlist) t)
     (declare (list nlist))
     (if (not (two-arg-<= n (car nlist))) (return nil))))

(defun >= (number &rest more-numbers)
  "Returns T if arguments are in strictly non-increasing order, NIL otherwise."
  (declare (optimize (safety 2))
	   (dynamic-extent more-numbers))
  (do* ((n number (car nlist))
	(nlist more-numbers (cdr nlist)))
       ((atom nlist) t)
     (declare (list nlist))
     (if (not (two-arg->= n (car nlist))) (return nil))))


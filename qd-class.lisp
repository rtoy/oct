(in-package "QD")

(defclass quad-double ()
  ((qd :initform +qd-zero+
       :reader qd-value
       :initarg :value
       :type %quad-double)))

(defmethod print-object ((qd quad-double) stream)
  (format stream "#q~/qdi::qd-format/" (qd-value qd)))

(defun qd-class-reader (stream subchar arg)
  (declare (ignore subchar arg))
  (make-instance 'quad-double :value (read-qd stream)))

;; Yow!  We redefine the #q reader that is in qd-io.lisp to read in
;; and make a real quad-double float, instead of the hackish
;; %quad-double.
(set-dispatch-macro-character #\# #\Q #'qd-class-reader)

(defmethod two-arg-+ ((a quad-double) (b quad-double))
  (make-instance 'quad-double :value (add-qd (qd-value a) (qd-value b))))

(defmethod two-arg-+ ((a quad-double) (b real))
  (make-instance 'quad-double :value (add-qd-d (qd-value a) (float b 1d0))))

(defmethod two-arg-+ ((a real) (b quad-double))
  (make-instance 'quad-double :value (add-d-qd (float a 1d0) (qd-value b))))

(defmethod two-arg-+ ((a real) (b real))
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

(defmethod two-arg-- ((a real) (b real))
  (cl:- a b))

(defmethod unary-minus ((a real))
  (cl:- a))

(defmethod unary-minus ((a quad-double))
  (make-instance 'quad-double :value (neg-qd (qd-value a))))

(defun - (number &rest more-numbers)
  (if more-numbers
      (do ((nlist more-numbers (cdr nlist))
	   (result number))
	  ((atom nlist) result)
         (declare (list nlist))
	 (setq result (- result (car nlist))))
      (unary-minus number)))


(defmethod two-arg-* ((a quad-double) (b quad-double))
  (make-instance 'quad-double :value (mul-qd (qd-value a) (qd-value b))))

(defmethod two-arg-* ((a quad-double) (b real))
  (make-instance 'quad-double :value (mul-qd-d (qd-value a) (float b 1d0))))

(defmethod two-arg-* ((a real) (b quad-double))
  (make-instance 'quad-double :value (mul-qd-d (qd-value b) (float a 1d0))))

(defmethod two-arg-* ((a real) (b real))
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

(defmethod two-arg-/ ((a real) (b real))
  (cl:/ a b))

(defmethod unary-divide ((a real))
  (cl:/ a))

(defmethod unary-divide ((a quad-double))
  (make-instance 'quad-double :value (div-qd +qd-one+ (qd-value a))))

(defun / (number &rest more-numbers)
  (if more-numbers
      (do ((nlist more-numbers (cdr nlist))
	   (result number))
	  ((atom nlist) result)
         (declare (list nlist))
	 (setq result (/ result (car nlist))))
      (unary-divide number)))

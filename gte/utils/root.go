package utils

// https://people.sc.fsu.edu/~jburkardt/c_src/fsolve/fsolve.html
// https://scipy.github.io/devdocs/reference/generated/scipy.optimize.fsolve.html

import (
	"errors"
	"fmt"
	"math"
)

const (
	p1    float64 = 0.1
	p25   float64 = 0.25
	p5    float64 = 0.5
	p05   float64 = 0.05
	p001  float64 = 0.001
	p0001 float64 = 0.0001

	DBL_EPSILON float64 = 2.220_446_049_250_313_08e-16
	DBL_MAX     float64 = 1.797_693_134_862_315_71e+308

	tolerance float64 = 1e-6
)

var (
	ErrImproper           = errors.New("improper input parameters")
	ErrNumberCallExceeded = errors.New("number of calls to fcn has reached or exceeded 200*(n+1)")
	ErrTolTooSmall        = errors.New("tol is too small. no further improvement in the approximate solution x is possible.")
	ErrBadProcess         = errors.New("iteration is not making good progress")
	ErrUnknown            = errors.New("error unknown")
)

// Function that takes at least one (possibly vector) argument, and returns a value of the same length.
type function func(x []float64) []float64

// returns the Euclidean norm of a vector.
func enorm(n int, x []float64) (value float64) {
	for i := 0; i < n; i++ {
		value += x[i] * x[i]
	}
	return math.Sqrt(value)
}

// computes the QR factorization of an M by N matrix.
func qrfac(m, n int, a []float64, lda int, pivot bool, ipvt []int, lipvt int, rdiag []float64, acnorm []float64) {
	var minmn, kmax int

	// compute the initial column norms and initialize sveral arrays
	wa := make([]float64, n)

	for j := 0; j < n; j++ {
		acnorm[j] = enorm(m, a[j*lda:])
		rdiag[j] = acnorm[j]
		wa[j] = rdiag[j]

		if pivot {
			ipvt[j] = j
		}
	}

	//reduce A to R with Householder transformations
	minmn = n
	if m < n {
		minmn = m
	}

	for j := 0; j < minmn; j++ {
		if pivot {
			//bring the column of largest norm into the pivot position
			kmax = j
			for k := j; k < n; k++ {
				if rdiag[kmax] < rdiag[k] {
					kmax = k
				}
			}

			if kmax != j {
				for i := 0; i < m; i++ {
					temp := a[j*lda]
					a[i+j*lda] = a[i+kmax*lda]
					a[i+kmax*lda] = temp
				}
				rdiag[kmax] = rdiag[j]
				wa[kmax] = wa[j]
				k := ipvt[j]
				ipvt[j] = ipvt[kmax]
				ipvt[kmax] = k
			}
		}

		// compute the Householder transformation to reduce the J-th column of A to a multiple of the J-th unit vector
		ajnorm := enorm(m-j, a[j+j*lda:])

		if ajnorm != 0 {
			if a[j+j*lda] < 0.0 {
				ajnorm *= -1
			}

			for i := j; i < m; i++ {
				a[i+j*lda] /= ajnorm
			}
			a[j+j*lda] += 1

			//apply the transformation to the remaining columns and update the norms
			for k := j + 1; k < n; k++ {
				sum := 0.0
				for i := j; i < m; i++ {
					sum += a[i+j*lda] * a[i+k*lda]
				}

				temp := sum / a[j+j*lda]
				for i := j; i < m; i++ {
					a[i+k*lda] = a[i+k*lda] - temp*a[i+j*lda]
				}

				if pivot && rdiag[k] != 0.0 {
					temp = a[j+k*lda] / rdiag[k]
					rdiag[k] *= math.Sqrt(max(0.0, 1.0-temp*temp))
					if p05*(rdiag[k]/wa[k])*(rdiag[k]/wa[k]) <= DBL_EPSILON {
						rdiag[k] = enorm(m-1-j, a[(j+1)+k*lda:])
						wa[k] = rdiag[k]
					}
				}
			}
		}

		rdiag[j] = -ajnorm
	}
}

// estimates an N by N Jacobian matrix using forward differences
func fdjac1(fn function, n int, x []float64, fvec []float64, fjac []float64, ldfjac, ml, mu int, epsfcn float64, wa1, wa2 []float64) {
	eps := math.Sqrt(max(epsfcn, DBL_EPSILON))
	msum := ml + mu + 1

	// computation of dense approximate jacobian
	if n <= msum {
		for j := 0; j < n; j++ {
			temp := x[j]
			h := eps * math.Abs(temp)

			if h == 0.0 {
				h = eps
			}

			x[j] = temp + h
			wa1 = fn(x)
			x[j] = temp

			for i := 0; i < n; i++ {
				fjac[i+j*ldfjac] = (wa1[i] - fvec[i]) / h
			}
		}
	} else {
		// computation of a banded approximate jacobian
		for k := 0; k < msum; k++ {
			for j := k; j < n; j += msum {
				wa2[j] = x[j]
				h := eps * math.Abs(wa2[j])
				if h == 0.0 {
					h = eps
				}
				x[j] = wa2[j] + h
			}

			wa1 = fn(x)

			for j := k; j < n; j += msum {
				x[j] = wa2[j]
				h := eps * math.Abs(wa2[j])
				if h == 0.0 {
					h = eps
				}

				for i := 0; i < n; i++ {
					if j-mu <= i && i <= j+ml {
						fjac[i+j*ldfjac] = (wa1[i] - fvec[i]) / h
					} else {
						fjac[i+j*ldfjac] = 0.0
					}
				}
			}
		}
	}
}

// combines Gauss-Newton and gradient for a minimizing step.
func dogleg(n int, r []float64, diag []float64, qtb []float64, delta float64, x []float64, wa1, wa2 []float64) {
	// calculate the Gauss-Newton direction
	jj := (n*(n+1))/2 + 1

	for k := 1; k <= n; k++ {
		j := n - k + 1
		jp1 := j + 1
		jj = jj - k
		l := jj + 1
		sum := 0.0

		for i := jp1; i <= n; i++ {
			sum = sum + r[l-1]*x[i-1]
			l = l + 1
		}

		temp := r[jj-1]
		if temp == 0.0 {

			l = j
			for i := 1; i <= j; i++ {
				temp = max(temp, math.Abs(r[l-1]))
				l = l + n - i
			}

			temp = DBL_EPSILON * temp
			if temp == 0.0 {
				temp = DBL_EPSILON
			}
		}
		x[j-1] = (qtb[j-1] - sum) / temp
	}

	// test whether the Gauss-Newton direction is acceptable
	for j := 0; j < n; j++ {
		wa1[j] = 0.0
		wa2[j] = diag[j] * x[j]
	}

	qnorm := enorm(n, wa2)
	if qnorm <= delta {
		return
	}

	// The Gauss-Newton direction is not acceptable
	// Calculate the scaled gradient direction
	l := 0
	valueR := 0.0
	for j := 0; j < n; j++ {
		temp := qtb[j]
		for i := j; i < n; i++ {
			if l == 0 {
				valueR = diag[6*n-1]
			} else {
				valueR = r[l-1]
			}
			if i == 0 {
				diag[3*n-1] += valueR * temp
			} else {
				wa1[i-1] += valueR * temp
			}

			l++
		}
		wa1[j] = wa1[j] / diag[j]
	}

	//calculate the norm of the scaled gradient and test for the special case in which the scaled gradient is zero
	gnorm := enorm(n, wa1)
	sgnorm := 0.0
	alpha := delta / qnorm

	// calculate the point along the scaled gradient at which the quadratic is minimized
	if gnorm != 0.0 {
		for j := 0; j < n; j++ {
			wa1[j] = (wa1[j] / gnorm) / diag[j]
		}

		l := 0
		for j := 0; j < n; j++ {
			sum := 0.0
			for i := j; i < n; i++ {
				sum = sum + r[l]*wa1[i]
				l = l + 1
			}
			wa2[j] = sum
		}
		temp := enorm(n, wa2)
		sgnorm = (gnorm / temp) / temp
		alpha = 0.0

		// if the scaled gradient direction os not acceptable
		// calculate the point along the dogleg at which the quadratic is minimized
		if sgnorm < delta {
			bnorm := enorm(n, qtb)
			temp = (bnorm / gnorm) * (bnorm / gnorm) * (sgnorm / delta)
			temp = temp - (delta/qnorm)*(sgnorm/delta)*(sgnorm/delta) + math.Sqrt(math.Pow(temp-(delta/qnorm), 2)+(1.0-(delta/qnorm)*(delta/qnorm))*(1.0-(sgnorm/delta)*(sgnorm/delta)))
			alpha = ((delta / qnorm) * (1.0 - (sgnorm/delta)*(sgnorm/delta))) / temp
		}
	}

	// form appropriate convex combination of the Gauss-Newton direction and the sclaed gradient direction
	temp := (1.0 - alpha) * min(sgnorm, delta)
	for j := 0; j < n; j++ {
		x[j] = temp*wa1[j] + alpha*x[j]
	}
}

// updates the Q factor after a rank one update of the matrix
func r1updt(m, n int, s []float64, u, v, w []float64) bool {
	var tau, tan, cs, sn, cotan float64

	// initialize the diagonal element pointer
	jj := (n*(2*m-n+1))/2 - (m - n)

	// move the nontrivial part of the last column of S into W
	l := jj
	for i := n; i <= m; i++ {
		w[i-1] = s[l-1]
		l = l + 1
	}

	//rotate the vector V into a multiple of the N-th unit vector in such a way that a spike is introduced into W
	nm1 := n - 1
	for j := n - 1; 1 <= j; j-- {
		jj = jj - (m - j + 1)
		w[j-1] = 0.0

		if v[j-1] != 0.0 {
			// determine a Givens rotation which eliminates the J-th element of V
			if math.Abs(v[n-1]) < math.Abs(v[j-1]) {
				cotan = v[n-1] / v[j-1]
				sn = p5 / math.Sqrt(p25+p25*cotan*cotan)
				cs = sn * cotan

				tau = 1.0
				if 1.0 < math.Abs(cs)*DBL_MAX {
					tau = 1.0 / cs
				}
			} else {
				tan = v[j-1] / v[n-1]
				cs = p5 / math.Sqrt(p25+p25*tan*tan)
				sn = cs * tan
				tau = sn
			}

			// apply the transformation to V and store the information necessary to recover the Givens rotatiuon
			v[n-1] = sn*v[j-1] + cs*v[n-1]
			v[j-1] = tau

			// apply the transformation to S and extend the spike in W
			l = jj
			for i := j; i <= m; i++ {
				temp := cs*s[l-1] - sn*w[i-1]
				w[i-1] = sn*s[l-1] + cs*w[i-1]
				s[l-1] = temp
				l = l + 1
			}
		}
	}

	// add the spike from the rank 1 update to W
	for i := 1; i <= m; i++ {
		w[i-1] = w[i-1] + v[n-1]*u[i-1]
	}

	// eliminate the spike
	sign := false

	for j := 1; j <= nm1; j++ {
		// determine a Givens rotation which eliminates the J-th element of the spike
		if w[j-1] != 0.0 {
			if math.Abs(s[jj-1]) < math.Abs(w[j-1]) {
				cotan = s[jj-1] / w[j-1]
				sn = p5 / math.Sqrt(p25+p25*cotan*cotan)
				cs = sn * cotan
				tau = 1.0
				if 1.0 < math.Abs(cs)*DBL_MAX {
					tau = 1.0 / cs
				}
			} else {
				tan = w[j-1] / s[jj-1]
				cs = p5 / math.Sqrt(p25+p25*tan*tan)
				sn = cs * tan
				tau = sn
			}

			// apply the transformation to s and reduce the spike in w
			l = jj
			for i := j; i <= m; i++ {
				temp := cs*s[l-1] + sn*w[i-1]
				w[i-1] = -sn*s[l-1] + cs*w[i-1]
				s[l-1] = temp
				tau = sn
			}

			w[j-1] = tau
		}

		// test for zero diagonal elements in the output s

		if s[jj-1] == 0.0 {
			sign = true
		}
		jj = jj + (m - j + 1)
	}

	// move W back into the last column of the output S
	l = jj
	for i := n; i <= m; i++ {
		s[l-1] = w[i-1]
		l = l + 1
	}

	if s[jj-1] == 0.0 {
		sign = true
	}
	return sign
}

// multiplies an M by N matrix A by the Q factor
func r1mpyq(m, n int, a []float64, lda int, v, w []float64) {
	var c, s float64

	// apply the first set of Givens rotations to A
	for j := n - 2; 0 <= j; j-- {
		if 1.0 < math.Abs(v[j]) {
			c = 1.0 / v[j]
			s = math.Sqrt(1.0 - c*c)
		} else {
			s = v[j]
			c = math.Sqrt(1.0 - s*s)
		}

		for i := 0; i < m; i++ {
			temp := c*a[i+j*lda] - s*a[i+(n-1)*lda]
			a[i+(n-1)*lda] = s*a[i+j*lda] + c*a[i+(n-1)*lda]
			a[i+j*lda] = temp
		}
	}

	// apply the second set of Givens rotations to A
	for j := 0; j < n-1; j++ {
		if 1.0 < math.Abs(w[j]) {
			c = 1.0 / w[j]
			s = math.Sqrt(1.0 - c*c)
		} else {
			s = w[j]
			c = math.Sqrt(1.0 - s*s)
		}

		for i := 0; i < m; i++ {
			temp := c*a[i+j*lda] + s*a[i+(n-1)*lda]
			a[i+(n-1)*lda] = -s*a[i+j*lda] + c*a[i+(n-1)*lda]
			a[i+j*lda] = temp
		}
	}
}

// constructs the standard form of Q from its factored form
func qform(m, n int, q []float64, ldq int) {
	var minmn int

	// zero out the upper triangle of Q in the first min(M,N) columns
	if m < n {
		minmn = m
	} else {
		minmn = n
	}

	for j := 1; j < minmn; j++ {
		for i := 0; i <= j-1; i++ {
			q[i+j*ldq] = 0.0
		}
	}

	// initialize remaning columns to those of the identity matrix
	for j := n; j < m; j++ {
		for i := 0; i < m; i++ {
			q[i+j*ldq] = 0.0
		}
		q[j+j*ldq] = 1.0
	}

	// accumulate Q from its factored form
	wa := make([]float64, m)

	for k := minmn - 1; 0 <= k; k-- {
		for i := k; i < m; i++ {
			wa[i] = q[i+k*ldq]
			q[i+k*ldq] = 0.0
		}
		q[k+k*ldq] = 1.0

		if wa[k] != 0 {
			for j := k; j < m; j++ {
				sum := 0.0
				for i := k; i < m; i++ {
					sum = sum + q[i+j*ldq]*wa[i]
				}
				temp := sum / wa[k]
				for i := k; i < m; i++ {
					q[i+j*ldq] = q[i+j*ldq] - temp*wa[i]
				}
			}
		}
	}
}

// finds a zero of a system of N nonlinear equations
func hybrd(fn function, x, diag, fjac, r, qtf, wa1, wa2, wa3, wa4 []float64) int {
	var msum int
	var iwa []int
	var delta, pnorm, fnorm1, actred, prered, ratio, xnorm float64

	n := len(x)
	ldfjac := n
	maxfev := 200 * (n + 1)
	code := 0
	ml := n - 1
	mu := n - 1
	epsfcn := 0.0
	factor := 100.0
	nfev := 0

	vec := make([]float64, n)

	// check the input parameters
	if n <= 0 || maxfev <= 0 || factor <= 0.0 || ldfjac < n {
		return code
	}

	// evaluate the function at the starting point and calculate ist norm
	vec = fn(x)
	nfev++
	fnorm := enorm(n, vec)

	// determine the number of calls to FCN needed to compute the jacobian matrix
	if ml+mu+1 < n {
		msum = ml + mu + 1
	} else {
		msum = n
	}

	// initialize iteration counter and monitors
	iter := 1
	ncsuc := 0
	ncfail := 0
	nslow1 := 0
	nslow2 := 0

	// beginning of the outer loop
	for {
		jeval := true

		// calculate the jacobian matrix
		fdjac1(fn, n, x, vec, fjac, ldfjac, ml, mu, epsfcn, wa1, wa2)

		nfev = nfev + msum

		// compute the QR factorization of the jacobian
		qrfac(n, n, fjac, ldfjac, false, iwa, 1, wa1, wa2)

		// on the first iteration and MODE is 1, scale according to the norms of the columns of the initial jacobian
		if iter == 1 {
			// on the first iteration, calculate the norm of the scaled X and initialize the step bound DELTA
			for j := 0; j < n; j++ {
				wa3[j] = diag[j] * x[j]
			}

			xnorm := enorm(n, wa3)
			if xnorm == 0.0 {
				delta = factor
			} else {
				delta = factor * xnorm
			}
		}

		// Form Q' * FVEC and store in QTF
		for i := 0; i < n; i++ {
			qtf[i] = vec[i]
		}

		for j := 0; j < n; j++ {
			if fjac[j+j*ldfjac] != 0.0 {
				sum := 0.0
				for i := j; i < n; i++ {
					sum = sum + fjac[i+j*ldfjac]*qtf[i]
				}
				temp := -sum / fjac[j+j*ldfjac]
				for i := j; i < n; i++ {
					qtf[i] = qtf[i] + fjac[i+j*ldfjac]*temp
				}
			}
		}

		// copy the triangular factor of the QR factorization into R
		// DO NOT ADJUST THIS LOOP, BECAUSE OF L
		for j := 1; j <= n; j++ {
			l := j
			for i := 1; i <= j-1; i++ {
				r[l-1] = fjac[(i-1)+(j-1)*ldfjac]
				l = l + n - i
			}

			r[l-1] = wa1[j-1]
			if wa1[j-1] == 0.0 {
				fmt.Println("hybrd: Matrix is singular")
			}
		}

		// Accumulate the orthogonal factor in FJAC
		qform(n, n, fjac, ldfjac)

		// beginning of the inner loop
	L:
		for {
			// determine the direction P
			dogleg(n, r, diag, qtf, delta, wa1, wa2, wa3)

			// Store the direction P and Z+P. Calculate the norm of P
			for j := 0; j < n; j++ {
				wa1[j] = -wa1[j]
				wa2[j] = x[j] + wa1[j]
				wa3[j] = diag[j] * wa1[j]
			}

			pnorm = enorm(n, wa3)

			// on the first iteration, adjust the initial step bound
			if iter == 1 {
				delta = min(delta, pnorm)
			}

			// evaluate the function at X+P and calculate its norm
			wa4 = fn(wa2)
			nfev = nfev + 1
			fnorm1 = enorm(n, wa4)

			// compute the scaled actual reduction
			if fnorm1 < fnorm {
				actred = 1.0 - (fnorm1/fnorm)*(fnorm1/fnorm)
			} else {
				actred = -1.0
			}

			// compute the scaled predicted reduction
			// DO NOT ADJUST THIS LOOP, BECAUSE OF L
			l := 1
			for i := 1; i <= n; i++ {
				sum := 0.0
				for j := i; j <= n; j++ {
					sum = sum + r[l-1]*wa1[j-1]
					l = l + 1
				}
				wa3[i-1] = qtf[i-1] + sum
			}

			temp := enorm(n, wa3)
			if temp < fnorm {
				prered = 1.0 - (temp/fnorm)*(temp/fnorm)
			} else {
				prered = 0.0
			}

			// compute the ratio of the actual to the predicted reduction
			if 0.0 < prered {
				ratio = actred / prered
			} else {
				ratio = 0
			}

			// update the step bound
			if ratio < p1 {
				ncsuc = 0
				ncfail = ncfail + 1
				delta = p5 * delta
			} else {
				ncfail = 0
				ncsuc = ncsuc + 1

				if p5 <= ratio || 1 < ncsuc {
					delta = max(delta, pnorm/p5)
				}

				if math.Abs(ratio-1.0) <= p1 {
					delta = pnorm / p5
				}
			}

			// on successful iteration, update X, FVEC and their norms
			if p0001 <= ratio {
				for j := 0; j < n; j++ {
					x[j] = wa2[j]
					wa2[j] = diag[j] * x[j]
					vec[j] = wa4[j]
				}
				xnorm = enorm(n, wa2)
				fnorm = fnorm1
				iter = iter + 1
			}

			// determin the progress of the iteration
			nslow1 = nslow1 + 1
			if p001 <= actred {
				nslow1 = 0
			}
			if jeval {
				nslow2 = nslow2 + 1
			}
			if p1 <= actred {
				nslow2 = 0
			}

			// test for convergence
			if delta <= tolerance*xnorm || fnorm == 0.0 {
				return 1
			}
			// test for termination and stringent tolerance
			if maxfev <= nfev {
				return 2
			}
			if p1*max(p1*delta, pnorm) <= DBL_EPSILON*xnorm {
				return 3
			}
			if nslow2 == 5 {
				return 4
			}
			// странное поведение scipy по сравнению с https://people.sc.fsu.edu/~jburkardt/c_src/fsolve/fsolve.html
			// руализация в scipy на FORTRAN: https://github.com/scipy/scipy/blob/v1.13.0/scipy/optimize/minpack/hybrd.f#L413
			// при nslow1==23 логика между этим решением и scipy сошлась, но это костыль!!!
			if nslow1 == 23 {
				return 5
			}

			// criterion for recalculating jacobian approximation by forward differences
			if ncfail == 2 {
				break L
			}

			// calculate the rank one modification to the jacobian and update QTF if necessary
			for j := 0; j < n; j++ {
				sum := 0.0
				for i := 0; i < n; i++ {
					sum = sum + fjac[i+j*ldfjac]*wa4[i]
				}

				wa2[j] = (sum - wa3[j]) / pnorm
				wa1[j] = diag[j] * ((diag[j] * wa1[j]) / pnorm)
				if p0001 <= ratio {
					qtf[j] = sum
				}
			}

			// compute the QR factorizatiom of the updated jacobian
			r1updt(n, n, r, wa1, wa2, wa3)
			r1mpyq(n, n, fjac, ldfjac, wa2, wa3)
			r1mpyq(1, n, qtf, 1, wa2, wa3)

			jeval = false
		}
	}
}

// Solver systems of nonlinear equations.
func Root(fn function, x []float64) ([]float64, error) {
	n := len(x)

	// check the input
	if n <= 0 {
		return nil, ErrImproper
	}

	lwa := (n * (3*n + 13)) / 2
	wa := make([]float64, lwa)
	for j := 0; j < n; j++ {
		wa[j] = 1.0
	}

	index := 6*n + (n*(n+1))/2

	code := hybrd(fn, x, wa, wa[index:], wa[6*n:], wa[n:], wa[2*n:], wa[3*n:], wa[4*n:], wa[5*n:])
	if err := getErrorFromCode(code); err != nil {
		return nil, err
	}

	return x, nil
}

func getErrorFromCode(code int) error {
	switch code {
	case 0:
		return ErrImproper
	case 1: // success
		return nil
	case 2:
		return ErrNumberCallExceeded
	case 3:
		return ErrTolTooSmall
	case 4, 5:
		return ErrBadProcess
	default:
		return ErrUnknown
	}
}

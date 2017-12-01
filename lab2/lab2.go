package lab2

import (
	"fmt"
	"math"

	"github.com/gonum/matrix/mat64"
	log "github.com/sirupsen/logrus"
)

var (
	a, b, a1, a2, a3, a4, b1, b2, b3, c1, c2, c3, d1, d2, d3, n1, n2, n3, k1, k2, p1, p2, q1, q2 float64
	n                                                                                            = 10
)

func SetEnvironment(_a, _b, _a1, _a2, _a3, _a4, _b1, _b2, _b3, _c1, _c2, _c3, _d1, _d2, _d3, _n1, _n2, _n3, _k1, _k2, _p1, _p2, _q1, _q2 float64) {
	a = _a
	b = _b
	a1 = _a1
	a2 = _a2
	a3 = _a3
	a4 = _a4
	b1 = _b1
	b2 = _b2
	b3 = _b3
	c1 = _c1
	c2 = _c2
	c3 = _c3
	d1 = _d1
	d2 = _d2
	d3 = _d3
	n1 = _n1
	n2 = _n2
	n3 = _n3
	k1 = _k1
	k2 = _k2
	p1 = _p1
	p2 = _p2
	q1 = _q1
	q2 = _q2
}

func SetN(_n int) {
	n = _n
}

func getSys() (k, p, q, f func(float64) float64, alpha, betta, gamma, delta, nu1, nu2 float64) {
	// -(ku')' + pu' + qu = f
	k = func(x float64) float64 {
		return b1*math.Pow(x, k1) + b2*math.Pow(x, k2) + b3
	}
	p = func(x float64) float64 {
		return c1*math.Pow(x, p1) + c2*math.Pow(x, p2) + c3
	}
	q = func(x float64) float64 {
		return d1*math.Pow(x, q1) + d2*math.Pow(x, q2) + d3
	}

	f = func(x float64) float64 {
		var res float64
		res += -k(x) * (a1*n1*(n1-1)*math.Pow(x, n1-2) + a2*n2*(n2-1)*math.Pow(x, n2-2) + a3*n3*(n3-1)*math.Pow(x, n3-2))
		res += -(b1*k1*math.Pow(x, k1-1) + b2*k2*math.Pow(x, k2-1)) * (a1*n1*math.Pow(x, n1-1) + a2*n2*math.Pow(x, n2-1) + a3*n3*math.Pow(x, n3-1))
		res += p(x) * (a1*n1*math.Pow(x, n1-1) + a2*n2*math.Pow(x, n2-1) + a3*n3*math.Pow(x, n3-1))
		res += q(x) * (a1*math.Pow(x, n1) + a2*math.Pow(x, n2) + a3*math.Pow(x, n3) + a4)
		return res
	}

	alpha = a1*math.Pow(a, n1) + a2*math.Pow(a, n2) + a3*math.Pow(a, n3) + a4
	betta = a1*n1*math.Pow(a, n1-1) + a2*n2*math.Pow(a, n2-1) + a3*n3*math.Pow(a, n3-1)
	gamma = a1*math.Pow(b, n1) + a2*math.Pow(b, n2) + a3*math.Pow(b, n3) + a4
	delta = -(a1*n1*math.Pow(b, n1-1) + a2*n2*math.Pow(b, n2-1) + a3*n3*math.Pow(b, n3-1))
	nu1 = 0
	nu2 = 0

	return
}

func tSys(kk, pp, qq, ff func(float64) float64, alpha, betta, gamma, delta, nnu1, nnu2 float64) (k, p, q, f, sigma1, sigma2, nu1, nu2 func(float64) float64) {
	// p0u'' + p1u' + p2u = f
	p0 := func(x float64) float64 { return -kk(x) }
	p1 := func(x float64) float64 { return pp(x) - der(kk, x) }
	p2 := qq

	r := func(x float64) float64 {
		return math.Exp(ig(func(x0 float64) float64 { return p1(x0) / p0(x0) }, a, x, 1000))
	}

	// (ku')' - qu = -f
	k = r
	p = func(float64) float64 { return 0 }
	q = func(x float64) float64 { return -p2(x) / p0(x) * r(x) }
	f = func(x float64) float64 { return -ff(x) / p0(x) * r(x) }

	sigma1 = func(x float64) float64 {
		return betta / alpha * k(x)
	}
	nu1 = func(x float64) float64 {
		return nnu1 / alpha * k(x)
	}

	sigma2 = func(x float64) float64 {
		return -delta / gamma * k(x)
	}
	nu2 = func(x float64) float64 {
		return nnu2 / gamma * k(x)
	}

	return
}

func Run() {
	var u = func(x float64) float64 {
		return a1*math.Pow(x, n1) + a2*math.Pow(x, n2) + a3*math.Pow(x, n3) + a4
	}
	var x []float64

	for i := 0; i <= n+1; i++ {
		x = append(x, a+float64(i)/float64(n+1)*(b-a))
	}

	k, p, q, f, alpha, betta, gamma, delta, nnu1, nnu2 := getSys()

	k, p, q, f, _, _, _, _ = tSys(k, p, q, f, alpha, betta, gamma, delta, nnu1, nnu2)
	for i := range x {
		x[i] = a + float64(i)/float64(n+1)*(b-a)
	}

	h := (b - a) / float64(n+1)
	aa := func(x float64) float64 {
		return 1 / h * (ig(k, x-h, x, 10) -
			ig(func(t float64) float64 { return q(t) * (x - t) * (t - (x - h)) }, x-h, x, 10))
	}

	dd := func(q func(float64) float64, x float64) float64 {
		return 1 / h / h * (ig(func(t float64) float64 { return (t - (x - h)) * q(t) }, x-h, x, 10) +
			ig(func(t float64) float64 { return ((x + h) - t) * q(t) }, x, x+h, 10))
	}

	var (
		A = mat64.NewDense(n+2, n+2, nil)
		c = mat64.NewDense(n+2, 1, nil)
		F = mat64.NewDense(n+2, 1, nil)
	)

	fmt.Println()
	for i := 1; i <= n; i++ {
		//progress visualisation
		//if i%(n/20) == 0 {
		fmt.Printf("\rCreating system. %6.2f%% Done", float64(i)/float64(n)*100)
		//}

		A.Set(i, i-1, aa(x[i]))
		A.Set(i, i, -aa(x[i+1])-aa(x[i])-h*h*dd(q, x[i]))
		A.Set(i, i+1, aa(x[i+1]))

		F.Set(i, 0, -h*h*dd(f, x[i]))
	}
	fmt.Println()

	A.Set(0, 0, -3*alpha-2*h*betta)
	A.Set(0, 1, 4*alpha)
	A.Set(0, 2, -alpha)
	F.Set(0, 0, 2*nnu1*h)

	A.Set(n+1, n-1, gamma)
	A.Set(n+1, n, -4*gamma)
	A.Set(n+1, n+1, 3*gamma+2*delta*h)
	F.Set(n+1, 0, 2*nnu2*h)

	// Show matrix
	//for i := 0; i <= n+1; i++ {
	//	for j := 0; j <= n+1; j++ {
	//		fmt.Printf("%7.2f ", A.At(i, j))
	//	}
	//	fmt.Printf("| %7.2f\n", F.At(i, 0))
	//}

	//fmt.Println(A)

	if err := c.Solve(A, F); err != nil {
		log.Error(err)
	}

	// nev'yazka
	//var d = &mat64.Dense{}
	//d.Mul(A, c)
	//d.Sub(d, F)
	//for i := 0; i <= n+1; i++ {
	//	fmt.Println(d.At(i, 0))
	//}

	// latex matrix
	fmt.Printf("\nn = %d\n\n", n)
	fmt.Println(`\begin{tabular}{| c | c | c | c | c |}
	\hline
	x & u & y & $\delta_y$ & $\varepsilon_y$ \\
	\hline`)
	for i, x := range x {
		if i%(n/10) != 0 {
			//continue
		}
		u, y := u(x), c.At(i, 0)
		fmt.Printf("\t%7.2f & %7.2f & %7.2f & %9.4f & %9.4f\\%% \\\\\n", x, u, y, math.Abs(u-y), math.Abs(u-y)/u*100)
	}
	fmt.Println("\t\\hline\n\\end{tabular}")
}

func der(f func(float64) float64, x float64) float64 {
	eps := 0.0001
	return (f(x+eps) - f(x-eps)) / (2 * eps)
}

func der2(f func(float64) float64, x float64) float64 {
	eps := 0.0001
	return (der(f, x+eps) - der(f, x-eps)) / (2 * eps)
}

func ig(f func(float64) float64, a, b float64, N int) float64 {
	N *= 2

	var I = f(a) + f(b)

	for i := 1; i < N; i++ {
		xi := a + float64(i)/float64(N)*(b-a)
		if i%2 == 1 {
			I += 4 * f(xi)
		} else {
			I += 2 * f(xi)
		}
	}

	I *= (b - a) / 3 / float64(N)

	return I
}

package lab1

import (
	"math"
	"strconv"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"

	"github.com/gonum/matrix/mat64"
	log "github.com/sirupsen/logrus"
)

var (
	a, b, a1, a2, a3, a4, b1, b2, b3, c1, c2, c3, d1, d2, d3, n1, n2, n3, k1, k2, p1, p2, q1, q2 float64
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

func Run() {
	k := func(x float64) float64 {
		return b1*math.Pow(x, k1) + b2*math.Pow(x, k2) + b3
	}
	p := func(x float64) float64 {
		return c1*math.Pow(x, p1) + c2*math.Pow(x, p2) + c3
	}
	q := func(x float64) float64 {
		return d1*math.Pow(x, q1) + d2*math.Pow(x, q2) + d3
	}

	alpha := a1*math.Pow(a, n1) + a2*math.Pow(a, n2) + a3*math.Pow(a, n3) + a4
	betta := a1*n1*math.Pow(a, n1-1) + a2*n2*math.Pow(a, n2-1) + a3*n3*math.Pow(a, n3-1)
	gamma := a1*math.Pow(b, n1) + a2*math.Pow(b, n2) + a3*math.Pow(b, n3) + a4
	delta := -(a1*n1*math.Pow(b, n1-1) + a2*n2*math.Pow(b, n2-1) + a3*n3*math.Pow(b, n3-1))

	f := func(x float64) float64 {
		var res float64
		res += -k(x) * (a1*n1*(n1-1)*math.Pow(x, n1-2) + a2*n2*(n2-1)*math.Pow(x, n2-2) + a3*n3*(n3-1)*math.Pow(x, n3-2))
		res += -(b1*k1*math.Pow(x, k1-1) + b2*k2*math.Pow(x, k2-1)) * (a1*n1*math.Pow(x, n1-1) + a2*n2*math.Pow(x, n2-1) + a3*n3*math.Pow(x, n3-1))
		res += p(x) * (a1*n1*math.Pow(x, n1-1) + a2*n2*math.Pow(x, n2-1) + a3*n3*math.Pow(x, n3-1))
		res += q(x) * (a1*math.Pow(x, n1) + a2*math.Pow(x, n2) + a3*math.Pow(x, n3) + a4)
		return res
	}

	var phi []func(float64) float64
	var x []float64

	n := 30
	for i := 0; i < n; i++ {
		var ii = i
		phi = append(phi, func(x float64) float64 {
			return math.Pow(x-a, float64(ii)) * math.Pow(x-b, 2)
		})
		x = append(x, a+float64(ii+1)/float64(n+1)*(b-a))
	}

	phi[0] = func(x float64) float64 {
		var A = b + (b-a)*gamma/((b-a)*delta+2*gamma)
		return (x - a) * (x - a) * (x - A)
	}

	phi[1] = func(x float64) float64 {
		var B = a + (a-b)*alpha/((b-a)*betta+2*alpha)
		return (x - b) * (x - b) * (x - B)
	}

	var u = func(x float64) float64 {
		return a1*math.Pow(x, n1) + a2*math.Pow(x, n2) + a3*math.Pow(x, n3) + a4
	}

	var (
		A  = mat64.NewDense(n, n, nil)
		c1 = mat64.NewDense(n, 1, nil)
		c2 = mat64.NewDense(n, 1, nil)
		F  = mat64.NewDense(n, 1, nil)
		//r = mat64.NewDense(n, 1, nil)
	)

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			A.Set(i, j, -k(x[i])*der2(phi[j], x[i])+(p(x[i])-der(k, x[i]))*der(phi[j], x[i])+q(x[i])*phi[j](x[i]))
		}
		F.Set(i, 0, f(x[i]))
	}

	if err := c1.Solve(A, F); err != nil {
		log.Error(err)
	}

	var un1 = func(x float64) float64 {
		var ans float64

		for i := 0; i < n; i++ {
			ans += phi[i](x) * c1.At(i, 0)
		}

		return ans
	}

	N := 10
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			A.Set(i, j, mul(func(x float64) float64 {
				return -k(x)*der2(phi[i], x) + (p(x)-der(k, x))*der(phi[i], x) + q(x)*phi[i](x)
			}, func(x float64) float64 {
				return -k(x)*der2(phi[j], x) + (p(x)-der(k, x))*der(phi[j], x) + q(x)*phi[j](x)
			}, a, b, N))
		}
		F.Set(i, 0, mul(func(x float64) float64 {
			return -k(x)*der2(phi[i], x) + (p(x)-der(k, x))*der(phi[i], x) + q(x)*phi[i](x)
		}, f, a, b, N))
	}

	if err := c2.Solve(A, F); err != nil {
		log.Error(err)
	}

	//r.Mul(A, c)
	//r.Sub(r, F)
	//var rr float64
	//for i := 0; i < n; i++ {
	//	rr += r.At(i, 0) * r.At(i, 0)
	//}
	//
	//rr = math.Sqrt(rr)
	//
	//log.Infof("r=%.4f", rr)

	var un2 = func(x float64) float64 {
		var ans float64

		for i := 0; i < n; i++ {
			ans += phi[i](x) * c2.At(i, 0)
		}

		return ans
	}

	var fn = func(x float64) float64 {
		return -k(x)*der2(u, x) + (p(x)-der(k, x))*der(u, x) + q(x)*u(x)
	}

	var t, ffn, ff []float64
	for i := 0; i < n; i++ {
		var x = a + float64(i+1)/float64(n+1)*(b-a)
		t = append(t, fn(x)-f(x))
		ffn = append(ffn, fn(x))
		ff = append(ff, f(x))
	}

	log.Info(t)
	log.Info(ff)
	log.Info(ffn)

	//fPol := plotter.NewFunction(f)
	//fPol.Color = color.RGBA{R: 255, G: 0, B: 0, A: 255}
	//fnPol := plotter.NewFunction(fn)
	//fnPol.Color = color.RGBA{R: 0, G: 255, B: 0, A: 255}

	pl, _ := plot.New()
	pl.X.Min, pl.X.Max = a, b

	plotutil.AddLinePoints(pl,
		"real", pointsFromFunc(u, 50),
		"kolok", pointsFromFunc(un1, 50),
		"mnmk", pointsFromFunc(un2, 50))

	//pl.Y.Min, pl.Y.Max = -2000, 2000

	//pl.Add(fPol)
	//pl.Add(fnPol)

	//pp := plotter.NewFunction(func(x float64) float64 {
	//	return x
	//})
	//pp.Color = color.RGBA{R: 0, G: 255, B: 0, A: 1}
	//

	log.Println("plotting")
	if err := pl.Save(5*vg.Inch, 5*vg.Inch, "plot"+strconv.FormatInt(int64(n), 10)+".png"); err != nil {
		log.Error(err)
	}
	log.Println("stopped plotting")
}

func pointsFromFunc(f func(float64) float64, n int) plotter.XYs {
	pts := make(plotter.XYs, n)
	for i := range pts {
		pts[i].X = a + float64(i+1)*(b-a)/float64(n+1)
		pts[i].Y = f(pts[i].X)
	}
	return pts
}

func der(f func(float64) float64, x float64) float64 {
	eps := 0.000001
	return (f(x+eps) - f(x-eps)) / (2 * eps)
}

func der2(f func(float64) float64, x float64) float64 {
	eps := 0.000001
	return (der(f, x+eps) - der(f, x-eps)) / (2 * eps)
}

func mul(f1, f2 func(float64) float64, a, b float64, N int) float64 {
	N *= 2

	var I = f1(a)*f2(a) + f1(b)*f2(b)

	for i := 1; i < N; i++ {
		xi := a + float64(i)/float64(N)*(b-a)
		if i%2 == 1 {
			I += 4 * f1(xi) * f2(xi)
		} else {
			I += 2 * f1(xi) * f2(xi)
		}
	}

	I *= (b - a) / 6 / float64(N)

	return I
}

package main

import (
	"fmt"

	"github.com/gonum/matrix/mat64"
)

func main() {
	// Physics params
	var (
		u1     float64 = 20. + 273.15
		u0     float64 = -20. + 273.15
		lambda float64 = 0.77
		cc     float64 = 830.
		ro     float64 = 1600.
		gamma  float64 = 7.
		R      float64 = 0.5
		T      float64 = 60 * 60
	)

	// Physics trans
	var (
		T1     = lambda / (cc * ro * R * R) * T
		gamma1 = R / lambda * gamma
	)

	// Numeric params
	var (
		N, M int     = 30, 30
		sig  float64 = 0.5
	)

	var tao, h float64 = T1 / float64(M), 1 / float64(N)

	var x []float64
	for i := 0; i <= N; i++ {
		x = append(x, float64(i)*h)
	}

	var dx = make([]float64, N+1)
	dx[0] = (x[1]*x[1] - x[0]*x[0]) / (2 * h)
	dx[N] = (x[N]*x[N] - x[N-1]*x[N-1]) / (2 * h)

	for i := 1; i < N; i++ {
		dx[i] = (x[i+1]*x[i+1] - x[i-1]*x[i-1]) / (4 * h)
	}

	var dp = make([]float64, N+1)
	//dp[0] = (x[1]*x[1] - x[0]*x[0]) / (2 * h)
	dp[N] = (x[N]*x[N] - x[N-1]*x[N-1]) / (2 * h)

	for i := 1; i < N; i++ {
		dp[i] = (x[i]*x[i] - x[i-1]*x[i-1]) / (2 * h)
	}

	var y = make([]float64, N+1)

	for i := 0; i <= N; i++ {
		y[i] = (u1 - u0) / u0
	}

	printArr(y)
	fmt.Println()

	for j := 0; j < M; j++ {
		b := make([]float64, N+2)
		c := make([]float64, N+2)
		d := make([]float64, N+2)
		phi := make([]float64, N+2)

		/* Simple */
		A := mat64.NewDense(N+1, N+1, nil)
		yy := mat64.NewDense(N+1, 1, nil)
		xx := mat64.NewDense(N+1, 1, nil)

		A.Set(0, 0, -3)
		A.Set(0, 1, 4)
		A.Set(0, 2, -1)
		yy.Set(0, 0, 0)

		for i := 1; i < N; i++ {
			// d[i] * v[i-1] + c[i] * v[i] + b[i] * v[i+1] = phi[i]
			d[i] = -sig * tao / (h * h) * x[i]
			b[i] = -sig * tao / (h * h) * x[i+1]
			c[i] = x[i] - b[i] - d[i]
			phi[i] = x[i]*y[i] + tao*(1-sig)/(h*h)*(x[i+1]*(y[i+1]-y[i])-x[i]*(y[i]-y[i-1]))

			A.Set(i, i-1, d[i])
			A.Set(i, i, c[i])
			A.Set(i, i+1, b[i])
			yy.Set(i, 0, phi[i])
		}

		A.Set(N, N, 3+2*h*gamma1)
		A.Set(N, N-1, -4)
		A.Set(N, N-2, 1)
		yy.Set(N, 0, 0)

		xx.Solve(A, yy)

		for i := range y {
			y[i] = xx.At(i, 0)
		}

		/* END Simple */

		/* Hard */
		//b[1] = sig * tao / (h * h) * dp[1]
		//c[1] = -dx[0]*0.5 - b[1]
		//phi[1] = -0.5*dx[0]*y[0] - (1-sig)*tao/(h*h)*dp[1]*(y[1]-y[0])
		//
		//for i := 1; i < N; i++ {
		//	// d[i] * v[i-1] + c[i] * v[i] + b[i] * v[i+1] = phi[i]
		//
		//	//d[i] = -sig * tao / (h * h) * x[i-1]
		//	//b[i] = -sig * tao / (h * h) * x[i+1-1]
		//	//c[i] = x[i-1] - b[i] - d[i]
		//	//phi[i] = x[i-1]*y[i-1] + tao*(1-sig)/(h*h)*(x[i+1-1]*(y[i+1-1]-y[i-1])-x[i-1]*(y[i-1]-y[i-1-1]))
		//
		//	d[i+1] = sig * tao / (h * h) * dp[i]
		//	b[i+1] = sig * tao / (h * h) * dp[i+1]
		//	c[i+1] = -dx[i] - (d[i+1] + b[i+1])
		//	phi[i+1] = -dx[i]*y[i] - tao*(1-sig)/(h*h)*(dp[i+1]*(y[i+1]-y[i])-dp[i]*(y[i]-y[i-1]))
		//}
		//
		//d[N+1] = sig * tao / (h * h) * dp[N]
		//c[N+1] = -sig*tao/h*gamma1*x[N] - 0.5*dx[N] - d[N+1]
		//phi[N+1] = (1-sig)*tao/h*gamma1*x[N]*y[N] - 0.5*dx[N]*y[N] + (1-sig)*tao/(h*h)*dp[N]*(y[N]-y[N-1])
		//
		//m := make([]float64, N+2)
		//w := make([]float64, N+2)
		//m[2] = -b[1] / c[1]
		//w[2] = phi[1] / c[1]
		//
		//for i := 2; i <= N; i++ {
		//	m[i+1] = -b[i] / (c[i] + d[i]*m[i])
		//	w[i+1] = (phi[i] - d[i]*w[i]) / (c[i] + d[i]*w[i])
		//}
		//
		//v := make([]float64, N+2)
		//
		//v[N+1] = (phi[N+1] - d[N+1]*w[N+1]) / (c[N+1] + d[N+1]*m[N+1])
		//
		//for i := N + 1; i > 1; i-- {
		//	v[i-1] = m[i]*v[i] + w[i]
		//}
		//
		//for i := 0; i <= N; i++ {
		//	y[i] = v[i+1]
		//}

		/* END Hard */

		printArr(y)
		fmt.Println()
	}
}

func printArr(arr []float64) {
	var u0 float64 = -20. + 273.15
	for _, x := range arr {
		fmt.Printf("%7.2f", (x+1)*u0-273.15)
	}
}

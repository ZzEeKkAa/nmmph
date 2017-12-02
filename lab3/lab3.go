package main

func main(){
	var T float64 = 60
	var N,M int = 10, 10
	var tao,h float64 = T/float64(M),1/float64(N)
	var sig float64 = 0.5

	var x []float64

	for i:=0; i<=N; i++{
		x = append(x,float64(i)*h)
	}

	var gamma float64 = 7.

	var dx =make([]float64,N+1)

	dx[0]= (x[1]*x[1]-x[0]*x[0])/(2*h)
	dx[N] = (x[N]*x[N]-x[N-1]*x[N-1])/(2*h)

	for i:=1;i<N; i++{
		dx[i] = (x[i+1]*x[i+1]-x[i-1]*x[i-1])/(4*h)
	}

	b:= make([]float64,N+1)
	c:= make([]float64,N+1)
	d:= make([]float64,N+1)

	b[1]=sig*tao/(h*h)*dx[1]
}
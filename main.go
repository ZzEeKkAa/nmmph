package main

import (
	"github.com/ZzEeKkAa/nmmph/lab2"
	"github.com/ZzEeKkAa/nmmph/lab1"
)

func main() {
	//lab1.SetEnvironment(0, 1, 1, 2, 4, 5, 2, 1, 1, 2, 1, 1, 1, 3, 1, 5, 3, 1, 0, 2, 0, 2, 1, 0)
	lab1.SetEnvironment(1, 4, 1, 2, 4, 5, 2, 1, 1, 2, 1, 1, 1, 3, 1, 5, 3, 1, 0, 2, 0, 2, 1, 0)
	lab1.Run()

	lab2.SetEnvironment(1, 4,
		1, 2, 4, 5,
			2, 1, 1,
				2, 1, 1,
					1, 3, 1,
						5, 3, 1,
							0, 2, 0, 2, 1, 0)

	//lab2.SetEnvironment(1, 4,
	//	1, 2, -3, 5,
	//		2, 1, 1,
	//			2, 1, 1,
	//				1, 3, 1,
	//					3, 2, 1,
	//						0, 2, 0, 2, 1, 0)
	//
	//lab2.SetN(100)
	//lab2.Run()
	//
	//lab2.SetN(100)
	//lab2.Run()
	//
	//lab2.SetN(1000)
	//lab2.Run()
}

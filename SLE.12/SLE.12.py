import numpy as np
A = []
b = []
Ae = [
        [3,1,1],
        [1,5,1],
        [1,1,7]
    ]

be = [
    [5],
    [7],
    [9]
]
A0 = [
        [0,2,3],
        [1,2,4],
        [4,5,6]
    ]

b0 =  [
    [13],
    [17],
    [32]
]
A1 =  [
        [15,1,1],
        [1,17,1],
        [1,1,19]
    ]

b1 =  [
    [17],
    [19],
    [21]
]
A2 =  [
        [-15,1,1],
        [1,-17,1],
        [1,1,-19]
    ]
b2 =  [
    [-17],
    [-19],
    [-21]
]
A3 =  [
        [-15,16,17],
        [18,-17,14],
        [17,18,-19]
    ]
b3 = [
    [17],
    [19],
    [21]
]
A4 =  [
        [15,14,14],
        [14,17,14],
        [14,14,19]
    ]

b4 =  [
    [17],
    [19],
    [21]
]
x1 = []
x = []
x1.append(np.linalg.solve(A0,b0))
x1.append(np.linalg.solve(A1,b1))
x1.append(np.linalg.solve(A2,b2))
x1.append(np.linalg.solve(A3,b3))
x1.append(np.linalg.solve(A4,b4))
with open("Matrix_Vector_results.txt", "a") as file:
    for i in range (5):
        file.write('3')
        file.write('\n')
        file.write(str(x1[i]).replace('[','').replace(']',''))
        file.write('\n')
with open("Matrix_Vector_results.txt", "a") as file:
    eps = [1e-3, 1e-6]
    for i in range(0,20,2):
       A_ = [[0 for t in range(i//2 + 4)] for _ in range((i)//2 + 4)]
       B_ = [0 for s  in range ((i//2) + 4)]
       for l in range (2):
       
            for j in range ((i)//2 + 4):
                if ( j != ((i)//2 - 1 + 4)):
                    B_[j] = -1
                else:
                    B_[j] = 1
                for k in range ((i)//2 + 4):
                    if (j == k):
                        A_[j][k] = 1 + eps[l]*13
                    elif (k < j):
                        A_[j][k] = 0 + eps[l] * 13
                    else:
                        A_[j][k] = -1 - eps[l]*13            
            A.append(A_)
            b.append(B_)
            x.append(np.linalg.solve(A[i+l],b[i+l]))
    k = 0
    for i in range (0,20,2):
       for j in range (2):
           file.write(str(k + 4))
           file.write('\n')
           file.write(str(x[i + j])[1:-1])
           file.write('\n')
       k = k + 1  
    
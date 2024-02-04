from termcolor import cprint

def print_matrix(xs, ys, M):
    '''Print a matrix object wit tagged rows and columns.
        xs -> list of row tags
        ys -> list of column tags
        M -> the matrix itself, made with numpy.
        Example:
        xs = ['a1','a1']
        xs = ['b1','b2']
            b1  b2
        a1  1   3
        a2  0   2
    '''
    
    s = 3

    print("        ", end="")
    for i in range(len(ys)):
        cprint(f"{ys[i]:3}", "red", end=" ")
    print()
    for i in range(len(xs) + 1):
        for j in range(len(ys) + 1):
            if j == 0:
                if i > 0:
                    cprint(xs[i-1], "red", end=" ")
                else:
                    print("  ", end="")
            if i == 0 or j == 0:
              cprint(f"{M[i][j]:3}", "blue", end=" ")
            else:
              print(f"{M[i][j]:3}", end=" ")
        print()

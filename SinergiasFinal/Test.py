
import numpy as np
import matplotlib.pyplot as plt
import time

def mult_martix(m1, m2):

    c = [[0 for col in range(len(m2[0]))] for row in range(len(m1))]

    for i in range(len(m1)):
        for j in range(len(m2[0])):
            for k in range(len(m2)):
                c[i][j] += m1[i][k] * m2[k][j]
                
    return c

def div_martix(m1, m2):

    c = [[0 for col in range(len(m2[0]))] for row in range(len(m1))]

    for i in range(len(m1)):
        for j in range(len(m1[0])):
            c[i][j] = m1[i][j] / m2[i][j]
                
    return c

def tran_matrix(m1):

    c = [[0 for col in range(len(m1))] for row in range(len(m1[0]))]

    for i in range(len(m1)):
        for j in range(len(m1[0])):
            c[j][i] = m1[i][j]

    return c

def dot_product(m1, m2):

    c = [[0 for col in range(len(m2[0]))] for row in range(len(m1))]

    for i in range(len(m1)):
        for j in range(len(m1[0])):
            c[i][j] = m1[i][j] * m2[i][j]

    return c

def add_matrix(m1,m2,add_mode):

    c = [[0 for col in range(len(m2[0]))] for row in range(len(m1))]

    for i in range(len(m1)):
        for j in range(len(m1[0])):
            c[i][j] = m1[i][j] + m2[i][j]*(-1+2*add_mode)

    return c

def thr_matrix(m1,threshold,cummulative):

    valid = False

    if(cummulative):
        c = 0
        for i in range(len(m1)):
            for j in range(len(m1[0])):
                c += abs(m1[i][j])
        if(c <= threshold):
            valid = True
    
    else:
        valid = True
        c = 0
        for i in range(len(m1)):
            for j in range(len(m1[0])):
                c = abs(m1[i][j])
                if(c > threshold):
                    valid = False
    
    return valid
        
in_time = time.time()
for n in range(max_iter):

    H_mult_num = mult_martix(tran_matrix(W),V)
    H_mult_den = mult_martix(mult_martix(tran_matrix(W),W),H)
    H_mult = div_martix(H_mult_num, H_mult_den)
    H = dot_product(H,H_mult)

    W_mult_num = mult_martix(V,tran_matrix(H))
    W_mult_den = mult_martix(mult_martix(W,H),tran_matrix(H))
    W_mult = div_martix(W_mult_num,W_mult_den)
    W = dot_product(W,W_mult)

    V_final = mult_martix(W,H)
    diff = add_matrix(V,V_final,False)

    if(thr_matrix(diff,threshold,cummulative)):
        print("GOOD w/", n+1 , "iterations.")
        break

out_time = time.time()
print("Time taken:", out_time-in_time, "seconds.")

f = [0]*len(V_final)

for chan in range(len(V_final)):
    f[chan] = plt.figure(chan)
    plt.plot(V[chan])
    plt.plot(V_final[chan])
    plt.plot(diff[chan])

plt.show()

# for chan in V_final:
#     plt.plot(chan)

# for chan in diff:
#     plt.plot(chan)

# plt.show()


# for r in V:
#    print(r, "INICIAL")

# for r in V_final:
#    print(r, "FINAL")

# for r in diff:
#    print(r, "DIFF")
import random
import matplotlib.pyplot as plt
from math import log, e, ceil
import numpy as np

n, h = 10, 2
ls, res = [], []
mx = 100


def uniform(max):
    return random.randint(1, max)


def expon(h):
    x = random.uniform(0, 1)
    return (-log(1 - x)) / h


def makevib(sz):
    vb = [[], [], [], [], []]
    for i in range(len(vb)):
        for j in range(sz):
            vb[i].append(expon(h))
        vb[i].sort()
    return vb


def makeVariationVb(sz):
    vb = []
    for j in range(sz):
        vb.append(uniform(h))
    vb.sort()
    return vb

# zxc = makeVariationVb(1000)
# for j in range(len(zxc)):
#     zxc[j] = round(zxc[j], 0)
# print(zxc.count(0))
# print(zxc.count(1))


def EmpUni(sizevb):
    tmp = makevib(sizevb)
    xmx = 0
    for i in range(5):
        for j in range(i + 1, 5):
            x = abs(max(tmp[i]) - max(tmp[j]))
            print(abs(max(tmp[i]) - max(tmp[j])))
            xmx = max(xmx, x)
    print(xmx / sizevb)

    vbmx = max(max(tmp[0]), max(tmp[1]), max(tmp[2]), max(tmp[3]), max(tmp[4])) + 1
    vbmn = min(min(tmp[0]), min(tmp[1]), min(tmp[2]), min(tmp[3]), min(tmp[4]))
    for z in range(5):
        empvb = tmp[z]
        empvb.sort()
        ys = []
        k = 0
        for i in empvb:
            if i <= empvb[0]:
                ys.append(0)
            elif i > empvb[0] and i <= empvb[len(empvb) - 1]:
                ys.append(k / len(empvb))
            else:
                ys.append(1)
            k += 1
        empvb.append(vbmx)
        ys.append(1)
        empvb.insert(0, vbmn)
        ys.insert(0, 0)
        plt.step(empvb, ys, label='EPDF' + str(z + 1))
    CDFX = []
    CDFY = []
    i = 0
    while i < vbmx:
        CDFX.append(i)
        CDFY.append(1 - e**(-h*i))
        i += 0.01
    plt.plot(CDFX, CDFY, label='PDF')
    plt.legend()
    plt.show()


def EmpDisc(sizevb):
    tmp = makevib(sizevb)
    xmx = 0
    for i in range(5):
        for j in range(i + 1, 5):
            x = abs(max(tmp[i]) - max(tmp[j]))
            print(abs(max(tmp[i]) - max(tmp[j])))
            xmx = max(xmx, x)
    print(xmx/sizevb)
    vbmx = max(max(tmp[0]), max(tmp[1]), max(tmp[2]), max(tmp[3]), max(tmp[4]))
    vbmn = min(max(tmp[0]), max(tmp[1]), max(tmp[2]), max(tmp[3]), max(tmp[4]))
    for z in range(5):
        empvb = tmp[z]
        empvb.sort()
        ys = []
        k = 0
        for i in empvb:
            if i <= empvb[0]:
                ys.append(0)
            elif empvb[0] < i <= empvb[len(empvb) - 1]:
                ys.append(k / len(empvb))
            else:
                ys.append(1)
            k += 1
        empvb.append(vbmx)
        ys.append(1)
        empvb.insert(0, vbmn)
        ys.insert(0, 0)
        plt.step(empvb, ys, label='ECDF' + str(z + 1))
    CDFX = [0]
    CDFY = [0]
    i = 0
    while i < vbmx:
        CDFX.append(i + 1)
        CDFY.append(i / (vbmx-1))
        i += 1
    plt.stem(CDFX, CDFY, label='CDF')
    plt.legend()
    plt.show()


def VariationRangeUni(sizevb):
    quan_ans = []
    tmp = makeVariationVb(sizevb)
    for quantile in [0.1, 0.5, 0.7]:
        k = ceil(quantile * (sizevb - 1))
        if k + 1 < quantile * sizevb:
            quan_ans.append(tmp[k+1])
        elif k + 1 == quantile * sizevb:
            quan_ans.append(tmp[k]*tmp[k+1]/2)
        elif k + 1 > quantile * sizevb:
            quan_ans.append(tmp[k])
        print('разность с теоретическим = ' + str(quan_ans[len(quan_ans)-1] - (-log(1-quantile)/h)))
    print(quan_ans)

def VariationRangeDisc(sizevb):
    quan_ans = []
    tmp = makeVariationVb(sizevb)
    for quantile in [0.1, 0.5, 0.7]:
        k = ceil(quantile * (sizevb - 1))
        if k + 1 < quantile * sizevb:
            quan_ans.append(tmp[k+1])
        elif k + 1 == quantile * sizevb:
            quan_ans.append(tmp[k]*tmp[k+1]/2)
        elif k + 1 > quantile * sizevb:
            quan_ans.append(tmp[k])
        print('разность с теоретическим = ' + str(quan_ans[len(quan_ans)-1] - (mx * quantile)))
    print(quan_ans)


def FreqPolyUni(sizevb):
    pdfmx = 0
    for i in range(5):
        tmp = makeVariationVb(sizevb)
        pdfmx = max(pdfmx, max(tmp))
        FreqX = tmp.copy()
        # for j in range (len(tmp)):
        #     tmp[j] = round(tmp[j], 0)
        # for j in range (len(FreqX)):
        #     FreqX[j] = round(FreqX[j], 0)
        FreqX = list(set(FreqX))
        FreqX.sort()
        FreqY = []
        for j in FreqX:
            FreqY.append(tmp.count(j) / len(tmp))
        plt.plot(FreqX, FreqY, label='EPMF' + str(i + 1))
    CDFX = []
    CDFY = []
    i = 0
    while i < pdfmx:
        CDFX.append(i)
        CDFY.append(h * (e ** (-h * i)))
        i += 0.01
    plt.plot(CDFX, CDFY, label='PMF')
    plt.legend()
    plt.show()

def NuFreqD(sizevb):
    pdfmx = 0
    conweigth = 100
    for i in range(5):
        tmp = makeVariationVb(sizevb)
        pdfmx = max(pdfmx, max(tmp))
        hist, bins = np.histogram(tmp, conweigth, weights=np.ones(sizevb)/sizevb/10)
        # plt.plot(bins[0:conweigth], hist, label= 'EPMF' + str(i + 1))
        plt.bar(bins[0:conweigth], hist, label='EHMF' + str(i + 1))
    CDFX = []
    CDFY = []
    i = 0
    while i < pdfmx:
        CDFX.append(i)
        CDFY.append(1 / conweigth/10)
        i += 1
    plt.plot(CDFX, CDFY, label='PMF')
    plt.legend()
    # plt.ylim((1 / conweigth / 10) - 0.001, (1 / conweigth / 10) + 0.001)
    plt.show()


def NuFreq(sizevb):
    pdfmx = 0
    for i in range(5):
        tmp = makeVariationVb(sizevb)
        pdfmx = max(pdfmx, max(tmp))
        hist, bins = np.histogram(tmp, 100, weights=np.ones(sizevb))
        # plt.plot(bins[0:100], hist, label= 'EPMF' + str(i + 1))
        plt.bar(bins[0:100], hist,  label='EHMF' + str(i + 1))
    CDFX = []
    CDFY = []
    i = 0
    while i < pdfmx:
        CDFX.append(i)
        CDFY.append(h * (e ** (-h * i)))
        i += 0.01
    # plt.plot(CDFX, CDFY, label='PMF')
    plt.legend()
    plt.show()


def FreqPolyDisc(sizevb):
    pdfmx = 0
    for i in range(5):
        tmp = makeVariationVb(sizevb)
        pdfmx = max(pdfmx, max(tmp))
        FreqX = tmp.copy()
        for j in range(len(tmp)):
            tmp[j] = round(tmp[j], 0)
        for j in range(len(FreqX)):
            FreqX[j] = round(FreqX[j], 0)
        FreqX = list(set(FreqX))
        FreqX.sort()
        FreqY = []
        for j in FreqX:
            FreqY.append(tmp.count(j) / len(tmp))
        plt.plot(FreqX, FreqY, label='EPMF' + str(i + 1))
    CDFX = []
    CDFY = []
    i = 0
    while i < pdfmx:
        CDFX.append(i)
        CDFY.append(1/pdfmx)
        i += 1
    plt.plot(CDFX, CDFY, label='PMF')
    plt.legend()
    plt.show()


def PrintFirstUni():
    xs, ys = [], []
    for i in range(n):
        tmp = expon(h)
        xs.append(tmp)
        ys.append(h * e ** (-h * tmp))
    plt.bar(xs,ys)
    plt.xlim(0,max(xs))
    plt.show()


def PrintFirstDisc():
    uniform(mx)
    for i in res:
        ls.append(1/mx)
    plt.bar(res, ls)
    plt.show()


def PrintSecondEMPUni():
    for i in [5, 10, 100, 1000, 10**5]:
        EmpUni(i)


def PrintSecondEMPDisc():
    for i in [5, 10, 100, 1000, 10**5]:
        EmpDisc(i)


def PrintSecondVarUni():
    for i in [5, 10, 100, 1000, 10**5]:
        VariationRangeUni(i)


def PrintSecondVarDisc():
    for i in [5, 10, 100, 1000, 10**5]:
        VariationRangeDisc(i)


def Factorial(n):
    res = 1
    for i in range(n):
        res *= i
    return res

# DZ 3


vb3 = makevib(100000)


def vb_moment(k, tmp):
    res_moment = [[], [], [], [], []]
    for i in range(0, 5):
        s = 0
        for j in range(0, len(tmp[i])):
            s += tmp[i][j]**k
        res_moment[i] = (s / len(tmp[i]))
    return res_moment


avg_vb = vb_moment(1, vb3)


def vb_mid_moment(k, avg, tmp):
    res_moment = [0, 0, 0, 0, 0]
    for i in range(0, 5):
        s = 0
        for j in range(0, len(tmp[i])):
            s += (tmp[i][j] - avg[i]) ** k
        res_moment[i] = (s / len(tmp[i]))
    return res_moment


zxczxc = [0.500003881068812, 0.4990051375326624, 0.5002954156795388, 0.49757915935785757, 0.500483155801488]
for i in zxczxc:
    print(1/i)

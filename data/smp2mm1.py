import os
from operator import itemgetter


def smp2mm1(filein, fileout):
    with open(filein, 'r') as fin, open(fileout, 'w') as fout:
        fin.readline()
        n = int(fin.readline().split()[0])
        nz = []
        has_diag = [False for i in range(n)]
        for line in fin.readlines():
            ele = line.split()
            if (ele[0] == '0'):
                break
            nz.append([int(ele[0]), int(ele[1]), ele[2]])
            if ele[0] == ele[1]:
                has_diag[int(ele[0]) - 1] = True
        for i in range(n):
            if not has_diag[i]:
                nz.append([i + 1, i + 1, '0'])
        print("% This MM1 formatted matrix file is converted by smp2mm0.py",
              file=fout)
        print('%10d%10d%10d' % (n, n, len(nz)), file=fout)
        for ele in sorted(nz, key=itemgetter(0, 1)):
            print('%10d%10d%30s' % (ele[0], ele[1], ele[2]), file=fout)


if __name__ == '__main__':
    if not os.path.exists('MM1'):
        os.mkdir("MM1")
    for f in sorted(os.listdir('SMP')):
        fi = 'SMP/' + f
        fo = 'MM1/' + f
        print('{} -> {}'.format(fi, fo))
        smp2mm1(fi, fo)

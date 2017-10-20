#!/usr/bin/env python
import sys

def parseItpFile(lista):
    parDict = {}
    for line in lista:
        line = line.split(';')[0]
        line = line.strip()
        if line and not line.startswith('#'):
            if line.startswith('[ '):
                flag = line.split()[1]
                if not parDict.has_key(flag):
                    parDict[flag] = []
            else:
                u = []
                t = line.split()
                for i in t:
                    try: v = eval(i)
                    except: v = i
                    u.append(v)
                parDict[flag].append(u)
    return parDict

def compareDicts(d1,d2):
    flags = [('pairs',2), ('bonds',2), ('angles',3), ('dihedrals',4)]

    for flag, l in flags:
        print "    ==> Comparing %s" % flag
        flagD1 = [x[:l+1] for x in d1[flag]]
        flagD2 = [x[:l+1] for x in d2[flag]]

        tempD1 = flagD1[:]
        tempD2 = flagD2[:]

        _tempD1 = tempD1[:]
        _tempD2 = tempD2[:]


        for item in tempD1:
            if str(item) in str(tempD2):
                try: _tempD1.remove(item)
                except: print "\tDuplicate item %s in %s itp2" % (`item`, flag)
                try: _tempD2.remove(item)
                except: print "\tDuplicate item %s in %s itp1" % (`item`, flag)
            else:
                t = item[:-1]
                t.reverse()
                itemRev = t + [item[-1]]
                if str(itemRev) in str(tempD2):
                    try: _tempD1.remove(item)
                    except: print "\tDuplicate item %s in %s itp2" % (`item`, flag)
                    try: _tempD2.remove(itemRev)
                    except: print "\tDuplicate item %s in %s itp1" % (`item`, flag)

        tD1 = _tempD1[:]
        tD2= _tempD2[:]

        if _tempD1 and _tempD2:
            # specially for dih since it may need to resort indexes
            i2l = {}
            for i2 in tD2:
                t = i2[:-1]
                t.sort()
                i2s = t + [i2[-1]]
                i2l[str(i2s)] = i2
            for i1 in tD1:
                t = i1[:-1]
                t.sort()
                i1s = t + [i1[-1]]
                if str(i1s) in str(i2l.keys()):
                    tD1.remove(i1)
                    tD2.remove(i2l[str(i1s)])

        if tD1:
            tD1.sort()
            print "    itp1: ", tD1
        if tD2:
            tD2.sort()
            print "    itp2: ", tD2

        if not tD1 and not tD2:
            print "        %s OK" % flag

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "ERROR: it needs 2 itp file as input"
        print "      compare_itps.py file1.itp file2.itp"
        sys.exit()
    itp1,itp2 = sys.argv[1:3]

    listItp1 = file(itp1).readlines()
    listItp2 = file(itp2).readlines()

    dictItp1 = parseItpFile(listItp1)
    dictItp2 = parseItpFile(listItp2)

    compareDicts(dictItp1, dictItp2)

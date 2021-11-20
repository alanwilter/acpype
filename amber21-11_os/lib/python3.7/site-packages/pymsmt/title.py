"""
This module was written for printing the title of each program in the package.
"""


def print_title(program, version, group='all'):
    totleng = 66
    shrtlen = 66 - 4

    print("*" * totleng)
    welcmstrg = "Welcome to use the " + program + " program"
    print("*", welcmstrg.center(shrtlen), "*")
    versnstrg = "Version " + version
    print("*", versnstrg.center(shrtlen), "*")
    print("*", "Author: Pengfei Li".center(shrtlen), "*")

    if group == 'merz':
        print("*", "Merz Research Group".center(shrtlen), "*")
        print("*", "Michigan State University".center(shrtlen), "*")
    elif group == 'shs':
        print("*", "Hammes-Schiffer Research Group".center(shrtlen), "*")
        print("*", "Yale University".center(shrtlen), "*")
    elif group == 'pl':
        print("*", "Li Research Group".center(shrtlen),"*")
        print("*", "Loyola University Chicago".center(shrtlen), "*")
    elif group == 'all':
        print("*", "Merz Research Group".center(shrtlen), "*")
        print("*", "Michigan State University".center(shrtlen), "*")
        print("*", "AND".center(shrtlen), "*")
        print("*", "Hammes-Schiffer Research Group".center(shrtlen), "*")
        print("*", "Yale University".center(shrtlen), "*")
        print("*", "AND".center(shrtlen), "*")
        print("*", "Li Research Group".center(shrtlen),"*")
        print("*", "Loyola University Chicago".center(shrtlen), "*")
    print("*" * totleng)


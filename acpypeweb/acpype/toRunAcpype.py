#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
import MySQLdb
from datetime import datetime
from datetime import datetime
acMaxFileSize = 1  # x Mb
acRefresh = 10  # seconds
acPrjDir = '/acpype'
mysqlAC = 'acpype_acpypejob'
toRunAcpype = 'toRunAcpype.py'
logAcpype = 'acpypeOut.log'
logOut = 'scriptsOut.log'
logError = 'scriptsError.log'
acpypeExe = '/usr/local/bin/acpype'
DATABASE_ENGINE = 'mysql'
DATABASE_NAME = 'acpype_acpypejob'
DATABASE_USER = 'lkagami'
DATABASE_PASSWORD = '3085luc9973'
DATABASE_HOST = ''
DATABASE_PORT = ''

def job_end(ferror, fout, switch_log, t0, logApp):
    print 'Total Running Time: %ss' % str(datetime.now() - t0)
    fix_chmod()
    if switch_log:
        ferror.close()
        fout.close()
    #cmd = 'cat %s %s %s > %s' % (logError,logOut,logAcpype,logLog)
    cmd = 'echo ----------%s---------- > %s; cat %s >> %s; echo ----------%s---------- >> %s; cat %s >> %s; echo ----------%s---------- >> %s; cat %s >> %s' % (logError, logLog, logError, logLog, logOut, logLog, logOut, logLog, logApp, logLog, logApp, logLog)
    os.system(cmd)

def sendEmail(jemail, status, db, mysqlTable, jobdir, label):
    if not jemail:
        os.system('touch no_email_defined')
        print 'no_email_defined'
        return
    from email.MIMEText import MIMEText #@UnresolvedImport
    cursor = MySQLdb.cursors.DictCursor(db)
    jobsql = '%' + os.path.basename(jobdir) + '%'
    query = "select * from %s where jobdir like '%s'" % (mysqlTable, jobsql)
    cursor_exec(cursor, query)
    cursorDict = cursor.fetchone()
    jname = cursorDict['jname']
    jdate = cursorDict['jdate']
    juser = cursorDict['juser']
    if status == 'Finished':
        msg = "Your job:\n\n  %s \n\nfrom\n\n  %s\n\nhas finished with success.\n\n%s" % (jname, jdate, label)
    if status == 'Cancelled':
        msg = "Your job:\n\n  %s \n\nfrom\n\n  %s\n\nwas cancelled.\n\n%s" % (jname, jdate, label)
    if status == 'Failed':
        msg = "Your job:\n\n  %s \n\nfrom\n\n  %s\n\nhas failed.\n\n%s" % (jname, jdate, label)
    msgMail = MIMEText(msg)
    msgMail['Subject'] = '%s Job Status' % mysqlTable.split('_')[0].upper()
    emailAdmin = '"%s" <%s>' % ADMINS[-1]
    msgMail['From'] = emailAdmin
    msgMail['To'] = jemail
    msgMail['Date'] = time.ctime()
    toAddresses = [jemail, ADMINS[-1][1]]
    #s = smtplib.SMTP('localhost')
    #s = smtplib.SMTP('mole.bio.cam.ac.uk', port = 587)
    s = smtplib.SMTP('ppsw.cam.ac.uk', port=25)
    try:
        #s.sendmail(emailAdmin, jemail, msgMail.as_string())
        s.sendmail(emailAdmin, toAddresses, msgMail.as_string())
        out = os.system('touch email_sent_to_%s' % juser)
        print "email_sent_to: %s - status %s" % (juser, out)
    except Exception, detail:
        out = os.system('touch email_sent_Failed_to_%s' % juser)
        print "email_sent_Failed_to: %s - status %s" % (juser, out)
        print "Reason for error:", detail
    s.close()


def job_end(ferror, fout, switch_log, t0, logApp):
    print 'Total Running Time: %ss' % str(datetime.now() - t0)
    fix_chmod()
    if switch_log:
        ferror.close()
        fout.close()
    #cmd = 'cat %s %s %s > %s' % (logError,logOut,logAcpype,logLog)
    cmd = 'echo ----------%s---------- > %s; cat %s >> %s; echo ----------%s---------- >> %s; cat %s >> %s; echo ----------%s---------- >> %s; cat %s >> %s' % (logError, logLog, logError, logLog, logOut, logLog, logOut, logLog, logApp, logLog, logApp, logLog)
    os.system(cmd)

paraDict = eval(sys.argv[1])

#{'args': u'-di dmp.pdb -c bcc -a gaff  ', 'charge_method': u'bcc', 'jdate': 'Tue Dec 30 19:08:34 2008', 'juser': u'alan', 'jobdir': u'/Users/alan/workspace/webapps/static/acpype/alan_20081230T190834', 'net_charge': None, 'file': 'dmp.pdb', 'multiplicity': None, 'jname': 'dmp', 'email': u'alanwilter@gmail.com', 'atom_type': u'gaff'}

switch_log = paraDict.get('switch_log', 0)
args = paraDict.get('args')
charge_method = paraDict.get('charge_method')
jdate = paraDict.get('jdate')
juser = paraDict.get('juser')
jobdir = paraDict.get('jobdir')
jname = paraDict.get('jname')
jemail = paraDict.get('email')

if switch_log:
    ferror = open(logError, 'a')
    sys.stderr = ferror
    fout = open (logOut, 'a')
    sys.stdout = fout
else:
    ferror = fout = None

print '_open_toRunAcpype_'

t0 = datetime.now()

print 'PID = ', os.getpid()

print '_parameters_'
print '\ncommand:', sys.argv[1:]
print '\nargs:', args

def getFilesAcpype():
    '''Create the output bundle results *.acpype'''
    acFolder = jname + '.acpype'
    cursor = MySQLdb.cursors.DictCursor(db)
    if os.path.exists(acFolder):
        cmd = 'zip -r %s_acpype.zip %s' % (jname, acFolder)
        os.system(cmd)
        out = os.system('touch acpype_zipped')
        print "acpype_zipped: status ", out
        jstatus = 'Finished'
        jobsql = '%' + os.path.basename(jobdir) + '%'
        query = "update %s set jstatus='%s' where jobdir like'%s'" % (mysqlAC, jstatus, jobsql)
        cursor_exec(cursor, query)
        os.system('touch get_files_DONE')
        print "get_files: status ", out
    else:
        os.system('touch get_files_FAILED')
        print "get_files: status ", 1
        job_failed(db, mysqlAC, jobdir, jemail, ferror, fout, switch_log, t0, logAcpype, 'ACPYPE Server')
    return

def run_acpype():
    if switch_log:
        cmd = 'nice -9 %s %s > %s 2>&1' % (acpypeExe, args, logAcpype)
    else:
        cmd = 'nice -9 %s %s' % (acpypeExe, args)
    out = os.system(cmd)
    os.system('touch run_acpype_DONE; touch DONE; cp /dev/null DONE')
    print "run_acpype: status ", out
    return out

## Main
if __name__ == '__main__':
    ## Run ACPYPE
    label = 'ACPYPE Server http://webapps.ccpn.ac.uk/acpype/status'
    out = run_acpype()
    if out == 0:
        getFilesAcpype()
        try:
            sendEmail(jemail, 'Finished', db, mysqlAC, jobdir, label)
        except:
            print 'ACPYPE OK: sendEmail Failed'
    else:
        job_failed(db, mysqlAC, jobdir, jemail, ferror, fout, switch_log, t0, logAcpype, label)
        sys.exit(1)

    job_end(ferror, fout, switch_log, t0, logAcpype)
    sys.exit(0)

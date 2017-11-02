#!/usr/bin/env python
# -*- coding: utf-8 -*-
acMaxFileSize = 1  # x Mb
acRefresh = 10  # seconds
acPrjDir = '/acpype'
mysqlAC = 'acpype_acpypejob'
toRunAcpype = 'toRunAcpype.py'
logAcpype = 'acpypeOut.log'
PYTHONCCPN = '/bin/python2.7'
DATABASE_USER = 'ccpndjango'
from django.shortcuts import render_to_response, render
from django.http import HttpResponseRedirect, HttpResponse
from django.contrib.auth import logout
from django.contrib.auth.decorators import login_required
from acpype.models import AcpypeJob, AcpypeJobForm
import os
import MySQLdb

def about(request):
    return render_to_response('about.html', {'user_': request.user, 'folder': 'acpype'})

def upload(request):
    if request.method == 'POST': # If the form has been submitted...
        juser = request.user
        form = AcpypeJobForm(request.POST, request.FILES) #, instance=job)
        if form.is_valid():
            #file = form.cleaned_data['file']
            parameters = form.cleaned_data
            parameters['juser'] = juser.username
            parameters['email'] = juser.email
            _results = handleUploadedAcpype(parameters) # parameters is a Dict
            #print _results
            return HttpResponseRedirect('status')
    else:
        form = AcpypeJobForm() # An unbound form

    return render_to_response('upload.html', {
        'user_': request.user, 'form': form, 'maxFileSize': acMaxFileSize,
        'folder': 'acpype'
    })

def logout_view(request):
    logout(request)
    return HttpResponseRedirect('acpype')

def status(request):
    runFlag = False
    user_obj = request.user
    _out = callStatus(user_obj.username, mysqlAC)
    #print _out
    if user_obj.username == DATABASE_USER:
        job_list = AcpypeJob.objects.all() #@UndefinedVariable
        admin = True
    else:
        job_list = AcpypeJob.objects.filter(juser=user_obj).exclude(jstatus='Deleted') #@UndefinedVariable
        admin = False
    if 'Running' in str(job_list) or 'Submitted' in str(job_list): runFlag = True
    return render_to_response('status.html',
            {'user_': user_obj, 'job_list': job_list, 'admin': admin,
             'runFlag': runFlag, 'refresh': acRefresh, 'folder': 'acpype'})

def call_status_func(request):
    user_obj = request.user
    if request.POST:
        func, out = callStatusFunc(request, mysqlAC)
        if func == 'View Log':
            #print out
            #return HttpResponse(out, mimetype="text/plain")
            return render_to_response('acpype/view_log.html',
            {'user_': user_obj, 'content': out.get('file'), 'jobId' : out.get('jobId'),
             'folder': 'acpype', 'link':out.get('link')})
        elif func == 'Download':
            response = HttpResponse(out.get('file'), mimetype="'application/zip'")
            response['Content-Disposition'] = 'attachment; filename=%s' % out.get('label')
            return response
        else:
            return HttpResponseRedirect('status')
    else:
        return HttpResponseRedirect('/acpype')

def handleUploadedAcpype(parameters):
    '''EXTERNAL: Upload a file ACPYPE and run it'''
    juser = parameters.get('juser')
    #jemail = parameters.get('email')
    file = parameters.get('file')
    fileType = file.content_type
    charge_method = parameters.get('charge_method')
    net_charge = parameters.get('net_charge')
    atom_type = parameters.get('atom_type')
    multiplicity = parameters.get('multiplicity')

    if (not os.path.isdir(acPrjDir)):
        os.mkdir(acPrjDir)
    os.chdir(acPrjDir)

    fileName = str(file.name)
    parameters['file'] = fileName

    jdate = datetime.datetime.now()
    ndate = jdate.ctime()
    parameters['jdate'] = ndate
    jname = os.path.splitext(fileName)[0].replace(' ', '_')
    parameters['jname'] = jname

    strNow = jdate.strftime("%Y%m%dT%H%M%S")
    userdir = juser + '_' + strNow
    jobdir = os.path.join(acPrjDir, userdir)
    parameters['jobdir'] = jobdir
    if (not os.path.isdir(jobdir)):
        os.mkdir(jobdir)

    destination = open(os.path.join(jobdir, fileName), 'wb+')
    for chunk in file.chunks():
        destination.write(chunk)
    destination.close()

    os.chdir(jobdir)

    if fileType == 'application/zip':
        os.mkdir(TEMPFOLDER)
        unzip(fileName, TEMPFOLDER)
    elif fileType == 'application/x-gzip':
        os.mkdir(TEMPFOLDER)
        untar(fileName, TEMPFOLDER)

    if net_charge != None: net_charge = '-n %i' % net_charge
    else: net_charge = ''

    if multiplicity != None: multiplicity = '-m %i' % multiplicity
    else: multiplicity = ''

    args = "-di %s -c %s -a %s %s %s" % (fileName, charge_method, atom_type,
                                         net_charge, multiplicity)
    parameters['args'] = args
    parameters['switch_log'] = nLog

    if RUN_METHOD == 'spawn':
        jpid = spawn(PYTHONCCPN, toRunAcpype, str(parameters))
    else:
        jline = '''jpid = spawn('%s', '%s', "%s")''' % \
            (PYTHONCCPN, toRunAcpype, str(parameters))
        tmpCron = os.path.join(jobdir, tmpCronFile)
        dict2Cron = {'tmpCron':tmpCron, 'jobdir':jobdir, 'mysqlServer':DATABASE_HOST,
                     'mysqlUser':DATABASE_USER, 'mysqlPass':DATABASE_PASSWORD, 'mysqlDb':DATABASE_NAME,
                     'mysqlTable':mysqlAC, 'jline':jline, 'jobsql': '%' + os.path.basename(jobdir) + '%'}
        writeCronScript(dict2Cron)
        addToCrontab(jobdir)
        jpid = 0

    jstatus = 'Submitted'
    db = MySQLdb.connect(host = DATABASE_HOST, user = DATABASE_USER, passwd = DATABASE_PASSWORD, db = DATABASE_NAME)
    cursor = MySQLdb.cursors.DictCursor(db)
    query = "insert into %s values ('%s','%s','%s','%s','%s','%s',%i)" \
               % (mysqlAC, juser, jname, fileName, ndate, jobdir, jstatus, jpid)
    cursor_exec(cursor, query)
    msg = "OK acpype %s started" % (jobdir)#, datetime.datetime.now(), RUN_METHOD)
    print msg
    logger.info(msg)
    return msg

def callStatus(juser, dbTable):
    '''EXTERNAL: Set and get job status'''
    ## Operate on mysql DB
    db = MySQLdb.connect(host = DATABASE_HOST, user = DATABASE_USER, passwd = DATABASE_PASSWORD, db = DATABASE_NAME)
    cursor = MySQLdb.cursors.DictCursor(db)
    if juser == DATABASE_USER:
        query = "select * from %s" % dbTable
    else:
        query = "select * from %s where juser='%s'" % (dbTable, juser)
    cursor_exec(cursor, query)
    jobs = cursor.fetchall()
    if len(jobs) > 0:
        for job in jobs:
            jpid, jobdir, jstatus = int(job['jpid']), job['jobdir'], job['jstatus']
            if RUN_METHOD == 'spawn':
                try: pid_test, dummy = os.waitpid(jpid, os.WNOHANG)
                except OSError: pid_test = 1
            else:
                pid_test = commands.getstatusoutput('ps %i' % jpid)[0]
            if pid_test == 0 and jstatus not in ['Deleted', 'Cancelled', 'Failed', 'Finished']:
                jstatus = 'Running'
                if dbTable == mysqlGrid:
                    file_aria_out = os.path.join(jobdir, logAria)
                    if os.path.exists(file_aria_out):
                        #jiter = jobStage(file_aria_out)
                        jiter = jobStage2(jobdir)
                        #query = "update %s set jstatus='%s', jiter='%i' where jpid='%i' and jobdir='%s'" % (dbTable,jstatus,jiter,jpid,jobdir)
                        query = "update %s set jstatus='%s', jiter='%i' where jobdir='%s'" % (dbTable, jstatus, jiter, jobdir)
                        cursor_exec(cursor, query)
                elif dbTable == mysqlISD:
                    jsample = '0'
                    #sampleLog = os.path.join(jobdir, isdSampleLogFile)
                    sampleLog = glob1(jobdir, '*_log')
                    #if os.path.exists(sampleLog):
                    if sampleLog:
                        jsample = nSamples(os.path.join(jobdir, sampleLog[0]))
                    #query = "update %s set jstatus='%s', jsamples='%s' where jpid='%i' and jobdir='%s'" % (dbTable,jstatus,jsample,jpid,jobdir)
                    query = "update %s set jstatus='%s', jsamples='%s' where jobdir='%s'" % (dbTable, jstatus, jsample, jobdir)
                    cursor_exec(cursor, query)
                else: #elif dbTable == mysqlAC:
                    #query = "update %s set jstatus='%s' where jpid='%i' and jobdir='%s'" % (dbTable,jstatus,jpid,jobdir)
                    query = "update %s set jstatus='%s' where jobdir='%s'" % (dbTable, jstatus, jobdir)
                    cursor_exec(cursor, query)
            elif jstatus in ['Submitted', 'Resumed'] and jpid == 0:
                pass
            elif (jstatus not in ['Deleted', 'Cancelled', 'Failed']):
                #print "AWSS 1", jstatus
                fname = os.path.join(jobdir, 'DONE') #'results.zip')
                if os.path.exists(fname):
                    jstatus = 'Finished'
                else:
                    print "AWSS 2", jstatus, fname
                    jstatus = 'Failed'
                #query = "update %s set jstatus='%s' where jpid='%i' and jobdir='%s'" % (dbTable,jstatus,jpid,jobdir)
                query = "update %s set jstatus='%s' where jobdir='%s'" % (dbTable, jstatus, jobdir)
                cursor_exec(cursor, query)
    return "OK callStatus at %s %s" % (datetime.datetime.now(), dbTable)

def callStatusFunc(request, dbTable): #(juser,jpid,jname,jdate,jobdir,func,REQUEST):
    '''EXTERNAL: Apply commands to jobs'''
    ## Operate on mysql DB
    db = MySQLdb.connect(host = DATABASE_HOST, user = DATABASE_USER, passwd = DATABASE_PASSWORD, db = DATABASE_NAME)
    cursor = MySQLdb.cursors.DictCursor(db)
    user = request.user.username
    jemail = request.user.email
    vars = request.POST
    juser = vars.get('juser')
    jdate = vars.get('jdate')
    jpid = int(vars.get('jpid', 0))
    jname = vars.get('jname')
    jobdir = vars.get('jobdir')
    func = vars.get('func')
    jobId = vars.get('jobId')
    addSamplesString = '+%i' % addSamples
    ## Execute appropriate function
    if (func == 'Cancel Job'):
        # for all
        cancelJob(dbTable, cursor, jobdir, jpid)
        return func, None
    elif (func == 'Delete'):
        # For all
        deleteJob(dbTable, cursor, jobdir, jpid)
        return func, None
    elif (func == 'Results'):
        # For CCPNGrid
        prefolder = jobdir.replace(WEBAPPSHOME, '')
        resultFiles = glob.glob('%s/Results_*.zip' % jobdir)
        resultFiles.sort(key = len)
        resultJnameReal = os.path.basename(resultFiles[0].replace('.zip', '')) # if mysql failed, jname will be 'test' on disk
        if resultJnameReal != 'Results_' + jname:
            print "WARN: failled to get jname: %s x %s" % (jname, resultJnameReal)
        resname = resultJnameReal
        results = {}

        if user == DATABASE_USER:
            jobId = '%s | %s | %s' % (juser, jname, jdate)
        else:
            jobId = '%s | %s' % (jname, jdate)

        results['inputAria'] = prefolder + '/' + 'inputAria.xml'
        results['jobId'] = jobId
        results['results_all'] = prefolder + '/' + resname + '.zip'
        results['results_pdbs'] = prefolder + '/' + resname + '_PDBS.zip'
        results['results_docs'] = prefolder + '/' + resname + '_DOCS.zip'
        results['results_figs'] = prefolder + '/' + resname + '_FIGS.zip'
        results['results_prj'] = prefolder + '/' + resname + '/Updated_CCPN_Prj_' + jname + '.zip'
        results['report_it'] = prefolder + '/' + resname + '/DOCS/last_it/report'
        results['report_wt'] = prefolder + '/' + resname + '/DOCS/refine/quality_checks'
        results['rama_it_all'] = prefolder + '/' + resname + '/FIGS/Ensemble_last_it_01.pdf'
        results['rama_ref_all'] = prefolder + '/' + resname + '/FIGS/Ensemble_refine_01.pdf'
        results['rama_it_res'] = prefolder + '/' + resname + '/FIGS/Ensemble_last_it_06.pdf'
        results['rama_ref_res'] = prefolder + '/' + resname + '/FIGS/Ensemble_refine_06.pdf'
        results['cns_analysis'] = prefolder + '/' + resname + '/CNS_Analysis.zip'
        pdbs = []
        pdbsw = []
        dir = jobdir + '/' + resname + '/PDBS/'
        for item in os.listdir(dir):
            if item.count('water.pdb'):
                pdbsw.append(dir.replace(WEBAPPSHOME, '') + item) #.replace('.pdb',''))
            else:
                pdbs.append(dir.replace(WEBAPPSHOME, '') + item) #.replace('.pdb',''))
        pdbs.sort(natcmp)
        pdbsw.sort(natcmp)
        results['pdbs'] = pdbs
        results['pdbsw'] = pdbsw
        return func, results
    elif (func == 'View Log'):
        # for all
        fname0 = os.path.join(jobdir, logLog)
        if os.path.exists(fname0):
            fname = fname0
        else:
            fname1 = fname0.replace(logLog, logAria)
            if os.path.exists(fname1):
                fname = fname1
            else:
                fname2 = fname1.replace(logAria, logOut)
                if os.path.exists(fname2):
                    fname = fname2
        link = fname.replace(WEBAPPSHOME, '')
        file = open(fname, 'rb')
        if user == DATABASE_USER:
            jobId = '%s | %s | %s' % (juser, jname, jdate)
            file = file.read()
        else:
            jobId = '%s | %s' % (jname, jdate)
            file = shrinkText(file.readlines())
        return func, {'file':file, 'jobId':jobId, 'link':link}
    elif (func == 'Try Rerun'):
        # for CCPNGrid
        os.umask(0)
        query = "select * from %s where jobdir='%s'" % (dbTable, jobdir)
        cursor_exec(cursor, query)
        round = int(cursor.fetchone()['round'] + 1)
        os.chdir(jobdir)
        if RUN_METHOD == 'spawn':
            jpid = spawn(PYTHONCCPN, TORUNARIA, 'round=%i' % round, 'switch_log=%i' % nLog, 'jemail=%s' % jemail, 'CONDOR=%i' % CONDOR)
        else:
            jline = '''jpid = spawn('%s','%s','round=%s','switch_log=%i','jemail=%s','CONDOR=%i')''' % \
                    (PYTHONCCPN, TORUNARIA, round, nLog, jemail, CONDOR)
            tmpCron = os.path.join(jobdir, tmpCronFile)
            dict2Cron = {'tmpCron':tmpCron, 'jobdir':jobdir, 'mysqlServer':DATABASE_HOST,
                         'mysqlUser':DATABASE_USER, 'mysqlPass':DATABASE_PASSWORD, 'mysqlDb':DATABASE_NAME,
                         'mysqlTable':mysqlGrid, 'jline':jline, 'jobsql': '%' + os.path.basename(jobdir) + '%'}
            writeCronScript(dict2Cron)
            addToCrontab(jobdir)
            jpid = 0
        jstatus = 'Resumed'
        query = "update %s set jstatus='%s', jpid='%i', round='%i' where jobdir='%s'" % (dbTable, jstatus, jpid, round, jobdir)
        cursor_exec(cursor, query)
        msg = "Job %s resumed" % (jobdir)#, datetime.datetime.now())
        print msg
        logger.info(msg)
        return func, None
    elif (func == 'Remove from DB'):
        # For all
        removeJobFromDB(dbTable, cursor, jobdir, jpid)
        return func, None
    elif (func == 'Download'):
        #Applies only for ACPYPE and ISD
        label = '%s_%s.zip' % (jname, os.path.basename(dictTableApp[dbTable]))
        fname = os.path.join(jobdir, label)
        file = open(fname, 'rb').read()
        msg = "Job %s downloaded" % (jobdir)#, datetime.datetime.now())
        print msg
        logger.info(msg)
        return func, {'file':file, 'label':label}
    elif (func == 'Report'):
        #Applies only for ISD
        fname = commands.getoutput('find %s/analysis -name "*_report.pdf"' % jobdir)
        label = os.path.basename(fname)
        file = open(fname, 'rb').read()
        link = fname.replace(WEBAPPSHOME, '')
        return func, {'file':file, 'label':label}
    elif (func == addSamplesString):
        os.umask(0)
        parameters = {}
        parameters['jobdir'] = jobdir
        parameters['jname'] = jname
        parameters['email'] = jemail
        parameters['switch_log'] = nLog
        parameters['addSamples'] = addSamples
        query = "select * from %s where jobdir='%s'" % (dbTable, jobdir)
        cursor_exec(cursor, query)
        #round = int( cursor.fetchone()['round'] + 1 )
        os.chdir(jobdir)
        if RUN_METHOD == 'spawn':
            jpid = spawn(PYTHONCCPN, toRunIsd, str(parameters))
        else:
            jline = '''jpid = spawn('%s', '%s', "%s")''' % \
                (PYTHONCCPN, toRunIsd, str(parameters))
            tmpCron = os.path.join(jobdir, tmpCronFile)
            dict2Cron = {'tmpCron':tmpCron, 'jobdir':jobdir, 'mysqlServer':DATABASE_HOST,
                         'mysqlUser':DATABASE_USER, 'mysqlPass':DATABASE_PASSWORD, 'mysqlDb':DATABASE_NAME,
                         'mysqlTable':mysqlISD, 'jline':jline, 'jobsql': '%' + os.path.basename(jobdir) + '%'}
            writeCronScript(dict2Cron)
            addToCrontab(jobdir)
            jpid = 0
        jstatus = 'Resumed'
        query = "update %s set jstatus='%s', jpid='%i' where jobdir='%s'" % (dbTable, jstatus, jpid, jobdir)
        cursor_exec(cursor, query)
        msg = "Job %s resumed" % (jobdir)#, datetime.datetime.now())
        print msg
        logger.info(msg)
        return func, None
    elif func.startswith('/'):
        # For CCPNGrid
        link = str(func)
        fn = os.path.abspath(WEBAPPSHOME + link)
        file = open(fn, 'r').read()
        return 'Report', {'file':file, 'jobId':jobId, 'link':link}
    # return results to view.py

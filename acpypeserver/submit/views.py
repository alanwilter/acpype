# -*- coding: utf-8 -*-
import os, os.path, pymysql.cursors, shutil
from django.shortcuts import render, render_to_response
from django.views.generic import CreateView, ListView
from django.urls import reverse_lazy
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from submit.models import Submission
from django.http import HttpResponse, HttpResponseRedirect
from django.contrib.auth.decorators import login_required
from django.contrib.auth import login, authenticate
from django.contrib.auth.forms import UserCreationForm
from django.shortcuts import render, redirect
from django import forms
from acpypeserver import settings as acpypesetting
from submit.models import Submission
from .tasks import process
from .forms import SignUpForm, SubmissionForm
from django.utils import timezone
from acpypeserver.celery import app
from django.contrib.auth.decorators import user_passes_test
from celery import uuid
from .serializers import SubmissionSerializer
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status as statushttp
from django.http import Http404
from rest_framework.decorators import api_view
from rest_framework.renderers import JSONRenderer

DATABASE_HOST = acpypesetting.DATABASES['default']['HOST']
DATABASE_USER = acpypesetting.DATABASES['default']['USER']
DATABASE_PASSWORD = acpypesetting.DATABASES['default']['PASSWORD']
DATABASE_NAME = acpypesetting.DATABASES['default']['NAME']


def home(request):
        return render(request, 'index.html')


class AuthRequiredMiddleware(object):

    def process_request(self, request):
        if not request.user.is_authenticated:
            return HttpResponseRedirect(reverse('/login/'))
            return None

@login_required
def Run(request):
    user_name = request.user.username
    if request.method == 'POST':
        form = SubmissionForm(request.POST, request.FILES)
        if form.is_valid():
            file = Submission(molecule_file=request.FILES['molecule_file'])
            file.juser = user_name
            file.jstatus = 'Running'
            file.save()
            molecule_file = request.FILES['molecule_file']
            cm = request.POST.get('charge_method')
            nc = request.POST.get('net_charge')
            ml = request.POST.get('multiplicity')
            at = request.POST.get('atom_type')
            mf = molecule_file.name
            mfs = str(mf)
            task_id = uuid()
            file.jcelery_id = task_id
            name = ((str(mfs)).split('_')[0])
            file.jname = ((str(name)).split('.')[0])
            file.save()
            process_task = process.apply_async((user_name, cm, nc, ml, at, mfs, task_id), task_id=task_id)
        else:
            return render(request, 'submit.html', locals())
    return HttpResponseRedirect('/status/')


def callStatusFunc(request):
    if request.method == 'POST':
        func = request.POST.get('func')
        jpid = request.POST.get('jpid')
        
        if func == 'download':
            db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
            cursor = pymysql.cursors.DictCursor(db)
            sql = "SELECT `jzipped` FROM `submit_submission` WHERE `jcelery_id`=%s"
            cursor.execute(sql, (jpid))
            jzipped = cursor.fetchone()
            db.close()
            zip_filename = jzipped['jzipped']
            zip_path = acpypesetting.MEDIA_ROOT
            os.chdir(acpypesetting.MEDIA_ROOT)
            zipfile = open(zip_filename, 'rb')
            response = HttpResponse(zipfile, content_type='application/zip')
            name_zipfile = ((str(zip_filename)).split('_')[3])
            response['Content-Disposition'] = 'attachment; filename={}'.format(name_zipfile)
            return response

        elif func == 'log':
            db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
            cursor = pymysql.cursors.DictCursor(db)
            sql = "SELECT `jlog`, `juser`, `molecule_file`, `date` FROM `submit_submission` WHERE `jcelery_id`=%s"
            cursor.execute(sql, (jpid))
            jlog = cursor.fetchone()
            fname = jlog['jlog']
            fuser = jlog['juser']
            fmol = jlog['molecule_file']
            fdata = jlog['date']
            fdata_str = fdata.strftime(' %a %b %d %H:%M %Y')
            os.chdir(acpypesetting.MEDIA_ROOT)
            pageFile = open(fname, "r")
            pageText = pageFile.read();
            pageFile.close()
            db.close()
            job = fuser + " | " + fmol + " | " + fdata_str
            return render_to_response('view_log.html', {'file':pageText, 'jobId':job})

        elif func == 'delete':
            db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
            cursor = pymysql.cursors.DictCursor(db)
            sql = "SELECT `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
            cursor.execute(sql, (jpid))
            user_folder = cursor.fetchone()
            folder_name = user_folder['usr_folder']
            job = Submission.objects.get(jcelery_id=jpid)
            job.jstatus = "Deleted"
            job.save()
            os.chdir(acpypesetting.MEDIA_ROOT)
            rmdir = str(acpypesetting.MEDIA_ROOT + "/" + folder_name)
            try:
                shutil.rmtree(rmdir)
            except:
                pass

        elif func == 'cancel':
            db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
            cursor = pymysql.cursors.DictCursor(db)
            sql = "SELECT `molecule_file` FROM `submit_submission` WHERE `jcelery_id`=%s"
            cursor.execute(sql, (jpid))
            molecule_file = cursor.fetchone()
            mol = molecule_file['molecule_file']
            app.control.revoke(jpid)
            job = Submission.objects.get(jcelery_id=jpid)
            job.jstatus = "Cancelled"
            job.save()
            os.chdir(acpypesetting.MEDIA_ROOT)

            if os.path.exists(mol):
                os.remove(mol)
            else:
                pass

        elif func == 'delete_db':
            db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
            cursor = pymysql.cursors.DictCursor(db)
            sql = "DELETE FROM `submit_submission` WHERE `jcelery_id`= %s"
            cursor.execute(sql, (jpid))
            db.commit()
            db.close()

    return HttpResponseRedirect('/status/')


class input(CreateView):

    template_name = 'submit.html'
    model = Submission
    fields = ('molecule_file', 'charge_method', 'net_charge', 'multiplicity', 'atom_type')


class status(ListView):
    model = Submission
    template_name = 'status.html'

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        data['is_running'] = Submission.objects.filter(juser=self.request.user, jstatus='Running').exists()
        return data

    def get_queryset(self):
        return Submission.objects.filter(juser=self.request.user).exclude(jstatus='Deleted').exclude(jstatus='Deleted_by_time')


class adminstatus(ListView):
    model = Submission
    template_name = 'status.html'

    def get_queryset(self):
        return Submission.objects.all()


def signup(request):
    if request.method == 'POST':
        form = SignUpForm(request.POST)
        if form.is_valid():
            form.save()
            username = form.cleaned_data.get('username')
            raw_password = form.cleaned_data.get('password1')
            user = authenticate(username=username, password=raw_password)
            login(request, user)
            return redirect('submit')
    else:
        form = SignUpForm()
    return render(request, 'signup.html', {'form': form})


"""
REST function

"""
@api_view(['GET', 'POST'])
def submit_list(request):
    
    if request.method == 'GET':

        user_name = request.GET.get('user_name')
        password = request.GET.get('password')
        func = request.GET.get('func')
        celery_id = request.GET.get('celery_id')
        mol_name = request.GET.get('molecule_file')
        user = authenticate(request, username=user_name, password=password)
        
        if user is not None:
            if func is None:
                submits = Submission.objects.filter(juser=user_name, molecule_file=mol_name)
                objects = SubmissionSerializer(submits, many=True)
                content = JSONRenderer().render(objects.data)
                return HttpResponse(content, content_type ='text/plain')

            elif func == 'gaff_mol2':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/'+mol_name+'_bcc_gaff.mol2'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    with open(fname) as f:
                        return HttpResponse(f, content_type ='text/plain')
                else:
                    return Response('Not found')

            elif func == 'pkl':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/'+mol_name+'.pkl'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    with open(fname) as f:
                        return HttpResponse(f, content_type ='text/plain')
                else:
                    return Response('Not found')

            elif func == 'ac_frcmod':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/'+mol_name+'_AC.frcmod'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    with open(fname) as f:
                        return HttpResponse(f, content_type ='text/plain')
                else:
                    return Response('Not found')

            elif func == 'ac_inpcrd':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/'+mol_name+'_AC.inpcrd'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    with open(fname) as f:
                        return HttpResponse(f, content_type ='text/plain')
                else:
                    return Response('Not found')

            elif func == 'ac_prmtop':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/'+mol_name+'_AC.prmtop'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    with open(fname) as f:
                        return HttpResponse(f, content_type ='text/plain')
                else:
                    return Response('Not found')

            elif func == 'charmm_inp':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/'+mol_name+'_CHARMM.inp'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    with open(fname) as f:
                        return HttpResponse(f, content_type ='text/plain')
                else:
                    return Response('Not found')

            elif func == 'ac_lib':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/'+mol_name+'_bcc_gaff.mol2'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    with open(fname) as f:
                        return HttpResponse(f, content_type ='text/plain')
                else:
                    return Response('Not found')

            elif func == 'charmm_prm':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/'+mol_name+'_CHARMM.prm'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    with open(fname) as f:
                        return HttpResponse(f, content_type ='text/plain')
                else:
                    return Response('Not found')

            elif func == 'charmm_rtf':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/'+mol_name+'_CHARMM.rtf'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    with open(fname) as f:
                        return HttpResponse(f, content_type ='text/plain')
                else:
                    return Response('Not found')

            elif func == 'cns_inp':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/'+mol_name+'_CNS.inp'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    with open(fname) as f:
                        return HttpResponse(f, content_type ='text/plain')
                else:
                    return Response('Not found')

            elif func == 'cns_par':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/'+mol_name+'_CNS.par'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    with open(fname) as f:
                        return HttpResponse(f, content_type ='text/plain')
                else:
                    return Response('Not found')

            elif func == 'cns_top':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/'+mol_name+'_CNS.top'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    with open(fname) as f:
                        return HttpResponse(f, content_type ='text/plain')
                else:
                    return Response('Not found')

            elif func == 'gmx_gro':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/'+mol_name+'_GMX.gro'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    readFile = open(fname, "r")
                    fileText = readFile.read();
                    readFile.close()
                    db.close()
                    print(fileText)
                else:
                    print('Not found')

            elif func == 'gmx_itp':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/'+mol_name+'_GMX.itp'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    with open(fname) as f:
                        return HttpResponse(f, content_type ='text/plain')
                else:
                    return Response('Not found')

            elif func == 'gmx_top':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/'+mol_name+'_GMX.top'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    with open(fname) as f:
                        return HttpResponse(f, content_type ='text/plain')
                else:
                    return Response('Not found')

            elif func == 'gmx_opls_itp':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/'+mol_name+'_GMX_OPLS.itp'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    with open(fname) as f:
                        return HttpResponse(f, content_type ='text/plain')
                else:
                    return Response('Not found')

            elif func == 'gmx_opls_top':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/'+mol_name+'_GMX_OPLS.top'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    with open(fname) as f:
                        return HttpResponse(f, content_type ='text/plain')
                else:
                    return Response('Not found')

            elif func == 'new_pdb':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/'+mol_name+'_NEW.pdb'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    with open(fname) as f:
                        return HttpResponse(f, content_type ='text/plain')
                else:
                    return Response('Not found')

            elif func == 'leap_log':
                db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
                cursor = pymysql.cursors.DictCursor(db)
                sql = "SELECT `molecule_file`, `usr_folder` FROM `submit_submission` WHERE `jcelery_id`=%s"
                cursor.execute(sql, (celery_id))
                jopen = cursor.fetchone()
                fmol = jopen['molecule_file']
                fold = jopen['usr_folder']
                mfs = str(fmol)
                name = ((str(mfs)).split('_')[0])
                mol_name = ((str(name)).split('.')[0])
                fname = fold+'/'+mol_name+'.acpype/leap.log'
                os.chdir(acpypesetting.MEDIA_ROOT)
                if os.path.exists(fname):
                    with open(fname) as f:
                        return HttpResponse(f, content_type ='text/plain')
                else:
                    return Response('Not found')
            else:
                return Response('Function not found')


        else:
            return Response(print(user_name),status=statushttp.HTTP_401_UNAUTHORIZED)

    elif request.method == 'POST':
        user_name = request.POST['user_name']
        password = request.POST['password']
        user = authenticate(request, username=user_name, password=password)
        if user is not None:
            form = SubmissionForm(request.POST, request.FILES)
            if form.is_valid():
                file = Submission(molecule_file=request.FILES['molecule_file'])
                file.juser = user_name
                file.jstatus = 'Running'
                file.save()
                molecule_file = request.FILES['molecule_file']
                cm = request.POST.get('charge_method')
                nc = request.POST.get('net_charge')
                ml = request.POST.get('multiplicity')
                at = request.POST.get('atom_type')
                mf = molecule_file.name
                mfs = str(mf)
                task_id = uuid()
                file.jcelery_id = task_id
                name = ((str(mfs)).split('_')[0])
                file.jname = ((str(name)).split('.')[0])
                file.save()
                process_task = process.apply_async((user_name, cm, nc, ml, at, mfs, task_id), task_id=task_id)
                return Response(status=statushttp.HTTP_202_ACCEPTED)
            else:
                Response(status=statushttp.HTTP_406_NOT_ACCEPTABLE)
        else:
           return Response(status=statushttp.HTTP_401_UNAUTHORIZED)
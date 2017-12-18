# -*- coding: utf-8 -*-
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
import os
from submit.models import Submission
from .tasks import process
import os.path
import pymysql.cursors
from .forms import SignUpForm, SubmissionForm
from django.utils import timezone
from acpypeserver.celery import app
from django.contrib.auth.decorators import user_passes_test

DATABASE_HOST = acpypesetting.DATABASES['default']['HOST']
DATABASE_USER = acpypesetting.DATABASES['default']['USER']
DATABASE_PASSWORD = acpypesetting.DATABASES['default']['PASSWORD']
DATABASE_NAME = acpypesetting.DATABASES['default']['NAME']


def home(request):
        return render(request,'index.html')

class AuthRequiredMiddleware(object):
    def process_request(self, request):
        if not request.user.is_authenticated:
            return HttpResponseRedirect(reverse('/login/'))
            return None

def Run(request):
    user_name = request.user.username
    if request.method == 'POST':
        form = SubmissionForm(request.POST, request.FILES)
        if form.is_valid():
            file = Submission(molecule_file = request.FILES['molecule_file'])
            file.juser = user_name
            file.jstatus = 'RUNNING'
            file.save()
            molecule_file = request.FILES['molecule_file']
            cm = request.POST.get('charge_method')
            nc = request.POST.get('net_charge')
            ml = request.POST.get('multiplicity')
            at = request.POST.get('atom_type')
            mf = molecule_file.name
            mfs = str(mf)
            process_task = process.delay(user_name,cm,nc,ml,at,mfs)
            file.jcelery_id = process_task.id
            name = ((str(mfs)).split('_')[0])
            file.jname = ((str(name)).split('.')[0])
            file.save()
        else:
            return render(request, 'submit.html', locals())
    return HttpResponseRedirect('/status/')

def callStatusFunc(request):
    if request.method == 'POST':
        func = request.POST.get('func')
        jpid = request.POST.get('jpid')
        db = pymysql.connect(host = DATABASE_HOST, user = DATABASE_USER, passwd = DATABASE_PASSWORD, db = DATABASE_NAME)
        cursor = pymysql.cursors.DictCursor(db)
        sql = "SELECT `jzipped` FROM `submit_Submission` WHERE `jcelery_id`=%s"
        cursor.execute(sql, (jpid))
        jzipped = cursor.fetchone()
        db.close()
        if func == 'download':
            zip_filename = jzipped['jzipped']+'.zip'
            zip_path = acpypesetting.MEDIA_ROOT
            os.chdir(acpypesetting.MEDIA_ROOT)
            zipfile = open(zip_filename, 'rb')
            response = HttpResponse(zipfile, content_type='application/zip')
            response['Content-Disposition'] = 'attachment; filename={}'.format(zip_filename)
            return response
        elif func == 'log':
            db = pymysql.connect(host = DATABASE_HOST, user = DATABASE_USER, passwd = DATABASE_PASSWORD, db = DATABASE_NAME)
            cursor = pymysql.cursors.DictCursor(db)
            sql = "SELECT `jlog` FROM `submit_Submission` WHERE `jcelery_id`=%s"
            cursor.execute(sql, (jpid))
            jlog = cursor.fetchone()
            fname = jlog['jlog']
            os.chdir(acpypesetting.MEDIA_ROOT)
            pageFile = open(fname, "r")
            pageText = pageFile.read();
            pageFile.close()
            db.close()
            return render_to_response('view_log.html', {'file':pageText, 'jobId':jpid})

        elif func == 'delete':
            job = Submission.objects.get(jcelery_id = jpid)
            job.jstatus = "DELETED"
            job.save()

        elif func == 'cancel':
            app.control.revoke(jpid)
            job = Submission.objects.get(jcelery_id = jpid)
            job.jstatus = "CANCELED"
            job.save()
            os.chdir(acpypesetting.MEDIA_ROOT)
            if os.path.exists(mfs):
                os.remove(mfs)
            else:
                pass
        elif func == 'delete_db':
            db = pymysql.connect(host = DATABASE_HOST, user = DATABASE_USER, passwd = DATABASE_PASSWORD, db = DATABASE_NAME)
            cursor = pymysql.cursors.DictCursor(db)
            sql = "DELETE FROM `submit_Submission` WHERE `jcelery_id`= %s"
            cursor.execute(sql, (jpid))
            db.commit()
            db.close()            
    return HttpResponseRedirect('/status/')

class input(CreateView):
    
    template_name = 'submit.html'
    model = Submission
    fields = ('molecule_file','charge_method','net_charge','multiplicity','atom_type')

class status(ListView):
    model = Submission
    template_name = 'status.html'
    def get_queryset(self):
        return Submission.objects.filter(juser=self.request.user).exclude(jstatus='Deleted')


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
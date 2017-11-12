# -*- coding: utf-8 -*-
from django.shortcuts import render, render_to_response
from django.views.generic import CreateView, ListView
from django.core.urlresolvers import reverse_lazy
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from submit.models import Submition
from django.http import HttpResponse, HttpResponseRedirect
from django.contrib.auth.decorators import login_required
from django.contrib.auth import login, authenticate
from django.contrib.auth.forms import UserCreationForm
from django.shortcuts import render, redirect
from django import forms
from acpypeserver import settings as acpypesetting
import os
from submit.models import Submition
from .tasks import process
from celery import shared_task
from celery.result import AsyncResult
import os.path
import MySQLdb
from .forms import SignUpForm
from django_celery_results.models import TaskResult
DATABASE_HOST = acpypesetting.DATABASES['default']['HOST']
DATABASE_USER = acpypesetting.DATABASES['default']['USER']
DATABASE_PASSWORD = acpypesetting.DATABASES['default']['PASSWORD']
DATABASE_NAME = acpypesetting.DATABASES['default']['NAME']


def home(request):
        return render(request,'index.html')

class SubmitionForm(forms.Form):

    molecule_file = forms.FileField(help_text='max. 1 megabyte', required=False)
    charge_method = forms.CharField(required=False)
    net_charge = forms.IntegerField(required=False)
    multiplicity = forms.IntegerField(required=False)
    atom_type = forms.CharField(required=False)

class AuthRequiredMiddleware(object):
    def process_request(self, request):
        if not request.user.is_authenticated():
            return HttpResponseRedirect(reverse('/login/'))
            return None

def submit_upload(request):
    if request.user.is_authenticated():
        return HttpResponseRedirect('/submit/')
    else:
        return HttpResponseRedirect('/login/')

def Run(request):
    if request.method == 'POST':
        form = SubmitionForm(request.POST, request.FILES)
        if form.is_valid():
            file = Submition(molecule_file = request.FILES['molecule_file'])
            file.save()
            molecule_file = request.FILES['molecule_file']
            cm = request.POST.get('charge_method')
            nc = request.POST.get('net_charge')
            ml = request.POST.get('multiplicity')
            at = request.POST.get('atom_type')
            mf = molecule_file.name
            mfs = str(mf)
            process.delay(cm,nc,ml,at,mfs)
            process_id = process.request.id
            tasks = AsyncResult(str(process_id)).state
            global tasks
        else:
            return render(request, 'submit.html', locals())
    return render(request, 'status.html')

def Download(request):
    
    db = MySQLdb.connect(host = DATABASE_HOST, user = DATABASE_USER, passwd = DATABASE_PASSWORD, db = DATABASE_NAME)
    cursor = MySQLdb.cursors.DictCursor(db)
    name1 = 'teste'
    zip_filename = name1+'_acpype.zip'
    zip_path = settings.MEDIA_ROOT
    os.chdir(settings.MEDIA_ROOT)
    zipfile = open(zip_filename, 'rb')
    response = HttpResponse(zipfile, content_type='application/zip')
    response['Content-Disposition'] = 'attachment; filename={}'.format(zip_filename)
    return response

class input(CreateView):
    
    template_name = 'submit.html'
    model = Submition
    fields = ('molecule_file','charge_method','net_charge','multiplicity','atom_type')


def call_status(request):
    mysqlAC = 'acpype_acpypejob'
    runFlag = False
    user_obj = request.user
    #_out = callStatus(user_obj.username, mysqlAC)
    #print _out
    if user_obj.username == DATABASE_USER:
        job_list = acpypeJob.objects.all() #@UndefinedVariable
        admin = True
    else:
        job_list = acpypeJob.objects.filter(juser=user_obj).exclude(jstatus='Deleted') #@UndefinedVariable
        admin = False
    if 'Running' in str(job_list) or 'Submitted' in str(job_list): runFlag = True
    return render_to_response('status.html',
            {'user_': user_obj, 'job_list': job_list, 'admin': admin,
             'runFlag': runFlag})

@login_required(login_url='/accounts/login/')
def status(request):
    template = "status.html"
    return render(request, template, {'tasks': tasks})
      

def signup(request):
    if request.method == 'POST':
        form = SignUpForm(request.POST)
        if form.is_valid():
            form.save()
            username = form.cleaned_data.get('username')
            raw_password = form.cleaned_data.get('password1')
            user = authenticate(username=username, password=raw_password)
            login(request, user)
            return redirect('home')
    else:
        form = SignUpForm()
    return render(request, 'signup.html', {'form': form})
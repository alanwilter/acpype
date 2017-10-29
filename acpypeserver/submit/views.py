# -*- coding: utf-8 -*-
from django.shortcuts import render, render_to_response
from django.views.generic import CreateView, ListView
from django.core.urlresolvers import reverse_lazy
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from submit.models import Submition
from django.http import HttpResponse
import os
from django import forms
from submit.models import Submition
import shutil
from acpypeserver import settings


def home(request):
        return render(request,'index.html')

class SubmitionForm(forms.Form):
    molecule_file = forms.FileField(help_text='max. 1 megabyte', required=False)
    charge_method = forms.CharField(required=False)
    net_charge = forms.IntegerField(required=False)
    multiplicity = forms.IntegerField(required=False)
    atom_type = forms.CharField(required=False)    

    def process(self):
        cm = self.cleaned_data.get('charge_method')
        nc = self.cleaned_data.get('net_charge')
        ml = self.cleaned_data.get('multiplicity')
        at = self.cleaned_data.get('atom_type')
        mf = self.cleaned_data.get('molecule_file')
        media_dir=os.chdir(settings.MEDIA_ROOT)
        execute_acpype = 'acpype -c {} -n {} -m {} -a {} -i {}'.format(cm,nc,ml,at,mf)
        os.system(execute_acpype)
        name = (str(mf)).split('.')[0]
        output_filename = '{}_acpype'.format(name)
        dir_name = '{}.acpype'.format(name)
        shutil.make_archive(output_filename, 'zip', dir_name)
        global name

    def Download(request):
        zip_filename = name+'_acpype.zip'
        zip_path = settings.MEDIA_ROOT
        zipfile = open(zip_filename, 'rb')
        response = HttpResponse(zipfile, content_type='application/zip')
        response['Content-Disposition'] = 'attachment; filename={}'.format(zip_filename)
        return response

class input(CreateView):
        template_name = 'submit.html'
        model = Submition
        fields = ('molecule_file','charge_method','net_charge','multiplicity','atom_type')

def Run(request):
    if request.method == 'POST':
        form = SubmitionForm(request.POST, request.FILES)
    if form.is_valid():
        file = Submition(molecule_file = request.FILES['molecule_file'])
        file.save()
        return render(request, 'display_results.html', using=form.process())
       
    else:
        return render(request, 'saved.html', locals())

def display_results(request):
    return render_to_response('display_results.html')




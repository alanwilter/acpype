# -*- coding: utf-8 -*-
from django.shortcuts import render, render_to_response
from django.views.generic import CreateView, ListView
from django.core.urlresolvers import reverse_lazy
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from submit.models import Submition
from submit.forms import SubmitionForm

def home(request):
        return render(request,'index.html')

class input(CreateView):
        template_name = 'submit.html'
        model = Submition
        fields = ('charge_method','net_charge','multiplicity','atom_type','molfile')

def Run(request):
    if request.method == 'POST':
        form = SubmitionForm(request.POST, request.FILES)
    if form.is_valid():
        file = Submition(molfile = request.FILES['molfile'])
        file.save()
        return render(request, 'display_results.html', using=form.process())
       
    else:
        return render(request, 'saved.html', locals())

def display_results(request):
    return render_to_response('display_results.html')

def Download(request):
    return render(request, 'index.html')


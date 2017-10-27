# -*- coding: utf-8 -*-
import os
from django import forms
from submit.models import Submition
from django.conf import settings
from django.http import HttpResponse
import shutil

class SubmitionForm(forms.Form):
    charge_method = forms.CharField(required=False)
    net_charge = forms.IntegerField(required=False)
    multiplicity = forms.IntegerField(required=False)
    atom_type = forms.CharField(required=False)
    molfile = forms.FileField(help_text='max. 2 megabytes', required=False)

    def process(self):
    	cm = self.cleaned_data.get('charge_method')
    	nc = self.cleaned_data.get('net_charge')
    	ml = self.cleaned_data.get('multiplicity')
    	at = self.cleaned_data.get('atom_type')
    	mf = self.cleaned_data.get('molfile')
    	media_dir=os.chdir("/home/server/teste2/acpypeserver/media")
    	execute_acpype = 'acpype -c {} -n {} -m {} -a {} -i {}'.format(cm,nc,ml,at,mf)
    	os.system(execute_acpype)
    	name = (str(mf)).split('.')[0]
    	output_filename = '{}_acpype'.format(name)
    	dir_name = '{}.acpype'.format(name)
    	shutil.make_archive(output_filename, 'zip', dir_name)


    

    




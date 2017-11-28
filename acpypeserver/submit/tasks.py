from __future__ import absolute_import, unicode_literals
import os
import shutil
from acpypeserver import settings
from celery import shared_task
import MySQLdb
import django.contrib.auth
from .models import Submition
from django.contrib.auth import login, authenticate
import os.path

@shared_task(ignore_result=False)
def process(user_name,cm,nc,ml,at,mfs):
	name_file = ((str(mfs)).split('.')[0])
	media_dir=os.chdir(settings.MEDIA_ROOT)
	os.chdir(settings.MEDIA_ROOT)
	execute_acpype = 'acpype -c {} -n {} -m {} -a {} -i {}'.format(cm,nc,ml,at,mfs)
	os.system(execute_acpype)
	output_filename = '{}_acpype'.format(name_file)
	dir_name = '{}.acpype'.format(name_file)
	shutil.make_archive(output_filename, 'zip', dir_name)
	job = Submition.objects.filter(juser=user_name).get(jstatus="STARTED")
	job.jstatus = "DONE"
	job.save()
	try:
		os.path.isfile(mfs)
		os.remove(mfs)
	except:
		pass
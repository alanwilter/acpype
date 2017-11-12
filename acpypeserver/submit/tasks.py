from __future__ import absolute_import, unicode_literals
import os
import shutil
from acpypeserver import settings
from celery import shared_task
import MySQLdb
import django.contrib.auth
DATABASE_HOST = settings.DATABASES['default']['HOST']
DATABASE_USER = settings.DATABASES['default']['USER']
DATABASE_PASSWORD = settings.DATABASES['default']['PASSWORD']
DATABASE_NAME = settings.DATABASES['default']['NAME']

@shared_task(ignore_result=False)
def process(cm,nc,ml,at,mfs):
	name_file = ((str(mfs)).split('.')[0])
	media_dir=os.chdir(settings.MEDIA_ROOT)
	os.chdir(settings.MEDIA_ROOT)
	execute_acpype = 'acpype -c {} -n {} -m {} -a {} -i {}'.format(cm,nc,ml,at,mfs)
	os.system(execute_acpype)
	output_filename = '{}_acpype'.format(name_file)
	dir_name = '{}.acpype'.format(name_file)
	shutil.make_archive(output_filename, 'zip', dir_name)
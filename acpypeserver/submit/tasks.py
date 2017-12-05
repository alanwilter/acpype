from __future__ import absolute_import, unicode_literals
import os
import shutil
from acpypeserver import settings
from celery import shared_task
import MySQLdb
import django.contrib.auth
from .models import Submition
import os.path
from time import gmtime, strftime
from django.core.mail import send_mail


DATABASE_HOST = settings.DATABASES['default']['HOST']
DATABASE_USER = settings.DATABASES['default']['USER']
DATABASE_PASSWORD = settings.DATABASES['default']['PASSWORD']
DATABASE_NAME = settings.DATABASES['default']['NAME']

@shared_task(ignore_result=False)
def process(user_name,cm,nc,ml,at,mfs):
	dt = strftime("%Y-%m-%d_%H:%M:%S", gmtime())
	dt_email = strftime("%Y-%m-%d %H:%M:%S", gmtime())
	name_file = ((str(mfs)).split('.')[0])
	media_dir=os.chdir(settings.MEDIA_ROOT)
	os.chdir(settings.MEDIA_ROOT)
	folder_name = user_name+'_'+dt
	execute_acpype = 'acpype -c {} -n {} -m {} -a {} -i {} -b {} > {}_{}.out'.format(cm,nc,ml,at,mfs,folder_name,user_name,dt)
	out = os.system(execute_acpype)
	if out == 0:
		output_filename = '{}_acpype_{}'.format(name_file,dt)
		dir_name = '{}.acpype'.format(folder_name)
		shutil.make_archive(output_filename, 'zip', dir_name)
		log_file = '{}_{}.out'.format(user_name,dt)
		job = Submition.objects.filter(juser=user_name).get(jstatus="STARTED")
		job.jstatus = "DONE"
		job.jzipped = output_filename
		job.jlog = log_file
		job.save()
		db = MySQLdb.connect(host = DATABASE_HOST, user = DATABASE_USER, passwd = DATABASE_PASSWORD, db = DATABASE_NAME)
		cursor = MySQLdb.cursors.DictCursor(db)
		sql = "SELECT `email` FROM `auth_user` WHERE `username`=%s"
		cursor.execute(sql, (user_name))
		get_email = cursor.fetchone()
		db.close()
		user_email = get_email['email']
		message = "Your Job '{}', has finished in {} \n\n ACPYPE Server Team ".format(name_file,dt_email)
		send_mail(
    	'ACPYPE Server',
    	message,
    	'from@example.com',
    	[user_email],
    	fail_silently=False,
		)
			
		try:
			os.path.isfile(mfs)
			os.remove(mfs)
		except:
			pass
	else:
		log_file = '{}_{}.out'.format(user_name,dt)
		job = Submition.objects.filter(juser=user_name).get(jstatus="STARTED")
		job.jstatus = "FAILED"
		job.jlog = log_file
		job.save()
		db = MySQLdb.connect(host = DATABASE_HOST, user = DATABASE_USER, passwd = DATABASE_PASSWORD, db = DATABASE_NAME)
		cursor = MySQLdb.cursors.DictCursor(db)
		sql = "SELECT `email` FROM `auth_user` WHERE `username`=%s"
		cursor.execute(sql, (user_name))
		get_email = cursor.fetchone()
		db.close()
		user_email = get_email['email']
		message = "Your Job '{}', has failed. \n\n ACPYPE Server Team ".format(name_file)
		send_mail(
    	'ACPYPE Server',
    	message,
    	'from@example.com',
    	[user_email],
    	fail_silently=False,
		)
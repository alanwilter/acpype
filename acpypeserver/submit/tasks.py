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
import yagmail


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
	os.system(execute_acpype)
	output_filename = '{}_acpype_{}'.format(name_file,dt)
	dir_name = '{}.acpype'.format(name_file)
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
	message = "Your Job '{}', has finished in {} \n ACPYPE Server Team ".format(name_file,dt_email)
	email = yagmail.SMTP('luciano8kagami@gmail.com', 'secret')
	email.send(user_email, 'ACPYPE Server', message)
			
	try:
		os.path.isfile(mfs)
		os.remove(mfs)
	except:
		pass
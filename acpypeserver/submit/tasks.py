from __future__ import absolute_import, unicode_literals
import shutil
from acpypeserver import settings
from celery import shared_task
import pymysql.cursors
import django.contrib.auth
from .models import Submission
import os.path
from datetime import datetime, timedelta
from django.core.mail import send_mail
from celery.task.schedules import crontab
from celery.decorators import periodic_task

DATABASE_HOST = settings.DATABASES['default']['HOST']
DATABASE_USER = settings.DATABASES['default']['USER']
DATABASE_PASSWORD = settings.DATABASES['default']['PASSWORD']
DATABASE_NAME = settings.DATABASES['default']['NAME']


@shared_task(ignore_result=False)
def process(user_name, cm, nc, ml, at, mfs):
    dt = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
    dt_email = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    name_file = ((str(mfs)).split('.')[0])
    media_dir = os.chdir(settings.MEDIA_ROOT)
    os.chdir(settings.MEDIA_ROOT)
    folder_name = user_name + '_' + dt
    execute_acpype = 'acpype -c {} -n {} -m {} -a {} -i {} -b {} > {}_{}.out'.format(cm, nc, ml, at, mfs, folder_name, user_name, dt)
    out = os.system(execute_acpype)
    if out == 0:
        output_filename = '{}_acpype_{}'.format(name_file, dt)
        dir_name = '{}.acpype'.format(folder_name)
        shutil.make_archive(output_filename, 'zip', dir_name)
        log_file = '{}_{}.out'.format(user_name, dt)
        job = Submission.objects.filter(juser=user_name).get(jstatus="RUNNING")
        job.jstatus = "FINISHED"
        job.jzipped = output_filename
        job.jlog = log_file
        job.save()
        try:
            shutil.rmtree(dir_name)
        except:
            pass
        if os.path.exists(mfs):
            os.remove(mfs)
        else:
            pass

        db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
        cursor = pymysql.cursors.DictCursor(db)
        sql = "SELECT `email` FROM `auth_user` WHERE `username`=%s"
        cursor.execute(sql, (user_name))
        get_email = cursor.fetchone()
        db.close()
        user_email = get_email['email']
        message = "Your Job '{}', has finished in {} \n\n ACPYPE Server Team ".format(name_file, dt_email)
        send_mail(
        'ACPYPE Server',
        message,
        'from@example.com',
        [user_email],
        fail_silently=False,
        )

    else:
        log_file = '{}_{}.out'.format(user_name, dt)
        job = Submission.objects.filter(juser=user_name).get(jstatus="RUNNING")
        job.jstatus = "FAILED"
        job.jlog = log_file
        job.save()
        try:
            shutil.rmtree(dir_name)
        except:
            pass
        if os.path.exists(mfs):
            os.remove(mfs)
        else:
            pass
        db = pymysql.connect(host=DATABASE_HOST, user=DATABASE_USER, passwd=DATABASE_PASSWORD, db=DATABASE_NAME)
        cursor = pymysql.cursors.DictCursor(db)
        sql = "SELECT `email` FROM `auth_user` WHERE `username`=%s"
        cursor.execute(sql, (str(user_name)))
        get_email = cursor.fetchone()
        db.close()
        user_email = get_email['email']
        message = "Your Job '{}', has failed. \n\n ACPYPE Server Team ".format(name_file)
        send_mail(
        'ACPYPE Server',
        message,
        'acpypeserver@gmail.com',
        [user_email],
        fail_silently=False,
        )


@periodic_task(run_every=(crontab(minute='*')), name="cleanup", ignore_result=True)
def cleanup():
    media_dir = os.chdir(settings.MEDIA_ROOT)
    os.chdir(settings.MEDIA_ROOT)
    init_date = datetime.today() - timedelta(days=14)
    final_date = datetime.today() - timedelta(days=7)
    job = Submission.objects.filter(date__range=[init_date, final_date]).values('jzipped')
    for file in job:
        if os.path.exists(file['jzipped'] + '.zip'):
            os.remove(file['jzipped'] + '.zip')
        else:
            pass
    job = Submission.objects.filter(date__range=[init_date, final_date]).values('jlog')
    for file in job:
        if os.path.exists(file['jlog']):
            os.remove(file['jlog'])
        else:
            pass
    try:
        for job in Submission.objects.filter(date__range=[init_date, final_date]):
            job.jstatus = 'DELETED'
            job.save()
    except:
        pass

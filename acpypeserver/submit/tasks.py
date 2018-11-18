from __future__ import absolute_import, unicode_literals
import shutil, os, os.path, csv, io, django.contrib.auth
from acpypeserver import settings
from celery import shared_task, app
from .models import Submission, MyUser
from datetime import datetime, timedelta
from django.core.mail import send_mail, EmailMessage
from celery.task.schedules import crontab
from celery.decorators import periodic_task
import re
import smtplib
import mimetypes
from email.mime.multipart import MIMEMultipart
from email import encoders
from email.message import Message
from email.mime.audio import MIMEAudio
from email.mime.base import MIMEBase
from email.mime.image import MIMEImage
from email.mime.text import MIMEText

DATABASE_HOST = settings.DATABASES['default']['HOST']
DATABASE_USER = settings.DATABASES['default']['USER']
DATABASE_PASSWORD = settings.DATABASES['default']['PASSWORD']
DATABASE_NAME = settings.DATABASES['default']['NAME']


@shared_task(ignore_result=False)
def process(user_name, cm, nc, ml, at, mfs, task_id):
    dt = datetime.now().strftime('%Y-%m-%d-%H:%M:%S')
    dt_email = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    name_file = ((str(mfs)).split('.')[0])
    media_dir = os.chdir(settings.MEDIA_ROOT)
    os.chdir(settings.MEDIA_ROOT)
    user_folder = user_name + '_' + dt + '_' + name_file
    os.makedirs(user_folder)
    path_to_molfile = settings.MEDIA_ROOT + '/' + mfs
    path_to_run = settings.MEDIA_ROOT + '/' + user_folder + '/'
    shutil.move(path_to_molfile, path_to_run)
    os.chdir(path_to_run)
    folder_name = name_file
    job = Submission.objects.filter(jcelery_id=task_id).get(jstatus="Queued")
    job.jstatus = "Running"
    job.charge_method = cm
    job.net_charge = nc
    job.multiplicity = ml
    job.atom_type = at
    job.save()
    execute_acpype = 'acpype -c {} -n {} -m {} -a {} -i {} -b {} > {}_{}.out'.format(cm, nc, ml, at, mfs, folder_name, user_name, dt)
    output_filename = '{}_acpype-{}'.format(name_file, dt)
    zip_name = output_filename + '.zip'
    log_file = '{}_{}.out'.format(user_name, dt)
    path_to_logfile = settings.MEDIA_ROOT + '/' + user_folder + '/' + log_file
    path_to_zipfile = settings.MEDIA_ROOT + '/' + user_folder + '/' + zip_name
    out = os.system(execute_acpype)

    if out == 0:
        dir_name = '{}.acpype'.format(folder_name)
        shutil.make_archive(output_filename, 'zip', dir_name)
        job.jstatus = "Finished"
        file = open (log_file, "r")
        lines = file.read()
        match = re.findall("Total time of execution:.*$",lines,re.MULTILINE)[0]
        spl = match.split(':')[1]
        job.runtime = spl
        job.jzipped = path_to_zipfile
        job.jlog = path_to_logfile
        job.usr_folder = user_folder
        job.save()
        eml = MyUser.objects.get(username=user_name)
        user_email = eml.email
        message = "Your Job '{}', has finished in {} \n\n ACPYPE Server Team ".format(name_file, dt_email)
        send_mail(
        'ACPYPE Server',
        message,
        'acpypeserver@gmail.com',
        [user_email],
        fail_silently=False,
        )

    else:
        log_file = '{}_{}.out'.format(user_name, dt)
        job = Submission.objects.filter(jcelery_id=task_id).get(jstatus="Running")
        job.jstatus = "Failed"
        job.jlog = path_to_logfile
        job.usr_folder = user_folder
        job.save()
        try:
            shutil.rmtree(dir_name)
        except:
            pass
        if os.path.exists(mfs):
            os.remove(mfs)
        else:
            pass
        eml = MyUser.objects.get(username=user_name)
        user_email = eml.email
        message = "Your Job '{}', has failed. \n\n ACPYPE Server Team ".format(name_file)
        send_mail(
        'ACPYPE Server',
        message,
        'acpypeserver@gmail.com',
        [user_email],
        fail_silently=False,
        )


@periodic_task(run_every=(crontab(minute=0, hour=0, day_of_week='sunday')), name="buildcsv", ignore_result=True)
def buildcsv():
    dt = datetime.now().strftime('%Y-%m-%d-%H:%M:%S')
    csv_name = dt + '_submit_table.csv'
    attachment_csv_file = io.StringIO()
    writer = csv.writer(attachment_csv_file)
    
    for obj in Submission.objects.all():
        writer.writerow([obj.jname, obj.charge_method,obj.net_charge,obj.multiplicity,obj.atom_type,obj.jstatus,obj.date,obj.runtime])

    csv_name2 = dt + '_user_table.csv'
    attachment_csv_file2 = io.StringIO()
    writer2 = csv.writer(attachment_csv_file2)
    
    for obj in MyUser.objects.all():
        writer2.writerow([obj.username, obj.first_name,obj.last_name,obj.email,obj.country])
    
      
    message = EmailMessage('Acpype Server - Tables', 'Tables', 'acpypeserver@gmail.com', ['acpypeserver@gmail.com'])
    message.attach(csv_name, attachment_csv_file.getvalue(), 'text/csv')
    message.attach(csv_name2, attachment_csv_file2.getvalue(), 'text/csv')
    message.send(fail_silently=False)

@periodic_task(run_every=(crontab(minute='*')), name="cleanup", ignore_result=True)
def cleanup():
    media_dir = os.chdir(settings.MEDIA_ROOT)
    os.chdir(settings.MEDIA_ROOT)
    init_date = datetime.today() - timedelta(days=14)
    final_date = datetime.today() - timedelta(days=7)
    job = Submission.objects.filter(date__range=[init_date, final_date]).values('usr_folder')

    for folder in job:

        if os.path.exists(folder['usr_folder']):
            shutil.rmtree(folder['usr_folder'])
        else:
            pass

    try:
        for job in Submission.objects.filter(date__range=[init_date, final_date]):
            job.jstatus = 'Deleted_by_time'
            job.save()
    except:
        pass

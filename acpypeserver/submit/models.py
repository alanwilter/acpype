from django.db import models
from django.utils import timezone

my_choices1 = [('bcc', "bcc (default)"), ('gas', "gasteiger"), ('user', "user")]
my_choices2 = [('gaff', "GAFF (default)"), ('gaff2', "GAFF2"), ('amber', "AMBER")]


class Submission(models.Model):
    jname = models.CharField(max_length=255, null=True, blank=True)
    molecule_file = models.FileField(upload_to='', null=True)
    charge_method = models.CharField(max_length=20, choices=my_choices1, blank=False, default='bcc')
    net_charge = models.IntegerField(default=0, null=True)
    multiplicity = models.PositiveIntegerField(default=1, null=True)
    atom_type = models.CharField(max_length=20, choices=my_choices2, blank=False, default='gaff')
    juser = models.CharField(max_length=30)
    jstatus = models.CharField(max_length=255, null=True, blank=True)
    jcelery_id = models.CharField(max_length=255, null=True, blank=True)
    jzipped = models.CharField(max_length=255, null=True, blank=True)
    jlog = models.CharField(max_length=255)
    date = models.DateTimeField(default=timezone.now)
    runFlag = models.CharField(max_length=255, default='False')
    usr_folder = models.CharField(max_length=255)

    class Meta:
        ordering = ['date']

        def __unicode__(self):
            return self.jcelery_id

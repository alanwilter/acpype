from django.db import models
from django.utils import timezone

my_choices1 = [('bcc', "bcc (default)"), ('gas', "gasteiger"), ('user', "user")]
my_choices2 = [('gaff', "GAFF (default)"), ('gaff2', "GAFF2"), ('amber', "AMBER")]


class Submission(models.Model):
    jname = models.CharField(max_length=255, null=True, blank=True)
    molecule_file = models.FileField(upload_to='', null=True, help_text="Required: Select a PDB, MDL or MOL2 file.")
    charge_method = models.CharField(max_length=20, choices=my_choices1, blank=False, default='bcc', help_text="Optional: Select a charge method.")
    net_charge = models.IntegerField(default=0, null=True, help_text="Optional: Enter an integer, if none ACPYPE will try to guess.")
    multiplicity = models.PositiveIntegerField(default=1, null=True, help_text="Optional: Integer (2S+1), default = 1.")
    atom_type = models.CharField(max_length=20, choices=my_choices2, blank=False, default='gaff', help_text="Optional: Select atom type.")
    juser = models.CharField(max_length=30)
    jstatus = models.CharField(max_length=255, null=True, blank=True)
    jcelery_id = models.CharField(max_length=255, null=True, blank=True)
    jzipped = models.CharField(max_length=255, null=True, blank=True)
    jlog = models.CharField(max_length=255,blank=True)
    date = models.DateTimeField(default=timezone.now)
    usr_folder = models.CharField(max_length=255, blank=True)

    class Meta:
        ordering = ['-date']

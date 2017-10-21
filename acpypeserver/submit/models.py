from django.db import models
from django.forms import ModelForm

my_choices1 = [(1, "bcc (default)"), (2, "gasteiger"), (5, "user")]
my_choices2 = [(1, "GAFF (default)"), (2, "AMBER")]

class Submition(models.Model):
        charge_method = models.CharField(choices=my_choices1, max_length=20, blank=False, default=1)
        net_charge = models.IntegerField(blank=False, default=1)
        multiplicity = models.IntegerField(blank=False, default=1)
        atom_type = models.CharField(choices=my_choices2, max_length=20, blank=False, default=1)
        
        
                      
class Upload(models.Model):
    description = models.CharField(max_length=255, blank=True)
    document = models.FileField(upload_to='home/')
    uploaded_at = models.DateTimeField(auto_now_add=True)
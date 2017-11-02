from django.db import models
from django import forms

acMaxFileSize = 1


class CommonJob(models.Model):
    juser = models.CharField(max_length=30)
    jname = models.CharField(max_length=255)
    fileName = models.CharField(max_length=255)
    jdate = models.CharField(max_length=30)
    jobdir = models.CharField(max_length=255, primary_key=True)
    jstatus = models.CharField(max_length=20)
    jpid = models.PositiveIntegerField()

    def __unicode__(self):
        return "%s %s %s %s" % (self.juser, self.jname, self.jdate, self.jstatus)

    class Meta:
        abstract = True
        ordering = ('juser', '-jobdir')  # , 'jname', 'jstatus','-jpid')


class AcpypeJob(CommonJob):
    pass


class AcpypeJobForm(forms.Form):
    file = forms.FileField(help_text="Required: Select a PDB, MDL or MOL2 file")
    #title = forms.CharField(max_length=255, required=False, help_text="(Optional)")
    charge_method = forms.ChoiceField(
        choices=(('bcc', 'bcc (default)'), ('gas', 'gasteiger'),
                 ('user', 'user')), required=False,
        help_text="Optional: Select a charge method")
    net_charge = forms.IntegerField(required=False, max_value=10, min_value=-10,
                                    help_text="Optional: enter an integer, if none ACPYPE will try to guess"
                                    )
    multiplicity = forms.IntegerField(required=False, max_value=10, min_value=-10,
                                      help_text="Optional: integer (2S+1), default = 1",
                                      initial='1')
    atom_type = forms.ChoiceField(
        choices=(('gaff', 'GAFF (default)'), ('amber', 'AMBER'),
                 #                                 ('bcc','BCC'), ('sybyl','SYBYL')
                 ), required=False,
        help_text='Optional: Select atom type'
    )

    def clean_file(self):
        file = self.cleaned_data['file']
        print file.content_type
        if len(file) > acMaxFileSize * 1024 * 1024:  # bytes
            raise forms.ValidationError('File size must not exceed %s Mb.' % acMaxFileSize)
        msg1 = 'File uploaded must be a valid PDB/MDL/MOL2 format.'
        if file.content_type in ['chemical/x-mol2', 'chemical/x-molfile', 'chemical/x-pdb']:
            return self.cleaned_data['file']
        elif file.content_type in ['application/octet-stream', 'text/plain'] and \
                file.name[-4:].upper() in ['MOL2', '.MDL', '.PDB', '.ENT']:
            return self.cleaned_data['file']
        else:
            raise forms.ValidationError(msg1)
        return self.cleaned_data['file']

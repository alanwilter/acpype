from django.contrib.auth.forms import AuthenticationForm 
from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User

acMaxFileSize = 1

class SignUpForm(UserCreationForm):
    first_name = forms.CharField(max_length=30, required=False, help_text='Optional.')
    last_name = forms.CharField(max_length=30, required=False, help_text='Optional.')
    email = forms.EmailField(max_length=254, help_text='Required. Inform a valid email address.')

    class Meta:
        model = User
        fields = ('username', 'first_name', 'last_name', 'email', 'password1', 'password2', )
        
class LoginForm(AuthenticationForm):
    username = forms.CharField(label="Username", max_length=30, 
                               widget=forms.TextInput(attrs={'class': 'form-control', 'name': 'username'}))
    password = forms.CharField(label="Password", max_length=30, 
                               widget=forms.PasswordInput(attrs={'class': 'form-control', 'name': 'password'}))

class SubmissionForm(forms.Form):
    molecule_file = forms.FileField()
    
    charge_method = forms.ChoiceField(
        choices=(('bcc', 'bcc (default)'), ('gas', 'gasteiger'),
                 ('user', 'user')), required=False,
        help_text="Optional: Select a charge method")
    net_charge = forms.IntegerField(required=False, max_value=10, min_value=-10)
    multiplicity = forms.IntegerField(required=False, max_value=10, min_value=-10, initial='1')
    atom_type = forms.ChoiceField(
        choices=(('gaff', 'GAFF (default)'), ('amber', 'AMBER'),
                 #                                 ('bcc','BCC'), ('sybyl','SYBYL')
                 ), required=False)
    runFlag = forms.CharField(max_length=255, required=False)
    usr_folder = forms.CharField(required=False)


    def clean_file(self):
        file = self.cleaned_data['molecule_file']
        print (file.content_type)
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
        return self.cleaned_data['molecule_file']
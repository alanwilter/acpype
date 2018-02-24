from django.contrib.auth.forms import AuthenticationForm
from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User
from submit.models import Submission
from django.core.validators import FileExtensionValidator
from acpypeserver import settings
from django.core.exceptions import ValidationError

def file_size(molecule_file): # add this to some file where you can import it from
        limit = settings.MAX_UPLOAD_SIZE * 1024 * 1024
        if molecule_file.size > limit:
            raise ValidationError('File too large. Size should not exceed 1 MiB.')

class SignUpForm(UserCreationForm):
    first_name = forms.CharField(max_length=30, required=False, help_text='Optional.')
    last_name = forms.CharField(max_length=30, required=False, help_text='Optional.')
    email = forms.EmailField(max_length=254, help_text='Required. Inform a valid email address.')

    class Meta:
        model = User
        fields = ('username', 'first_name', 'last_name', 'email', 'password1', 'password2',)


class LoginForm(AuthenticationForm):
    username = forms.CharField(label="Username", max_length=30,
                               widget=forms.TextInput(attrs={'class': 'form-control', 'name': 'username'}))
    password = forms.CharField(label="Password", max_length=30,
                               widget=forms.PasswordInput(attrs={'class': 'form-control', 'name': 'password'}))


class SubmissionForm(forms.Form):


    molecule_file = forms.FileField(validators=[FileExtensionValidator(settings.CONTENT_TYPES), file_size])

    charge_method = forms.ChoiceField(
            choices=(("bcc", 'bcc (default)'), ("gas", 'gasteiger'),
                     ("user", 'user')), required=False, help_text="Optional: Select a charge method.")
    net_charge = forms.IntegerField(required=False, max_value=10, min_value=-10, help_text="Optional: enter an integer, if none ACPYPE will try to guess.")
    multiplicity = forms.IntegerField(required=False, max_value=10, min_value=-10, initial='1', help_text="Optional: integer (2S+1), default = 1.")
    atom_type = forms.ChoiceField(
            choices=(("gaff", 'GAFF (default)'), ("gaff2", 'GAFF2'), ("amber", 'AMBER')),
            required=False,  help_text='Optional: Select atom type')
    runFlag = forms.BooleanField(required=False)
    usr_folder = forms.CharField(required=False)

    class Meta:
        model = Submission
        
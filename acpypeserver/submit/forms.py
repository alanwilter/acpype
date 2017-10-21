# -*- coding: utf-8 -*-
from django import forms
from submit.models import Submition
from submit.models import Upload

class SubmitionForm(forms.Form):

        class Meta:
                model = Submition


class UploadForm(forms.Form):
    class Meta:
        model = Upload
        fields = ('description', 'document', )

